# =============================================================================
# 05_ml_classifier.R  (Role 3 owns this)
# Frozen model handed to Stage 6/7 is trained on ALL discovery data (standard
# practice for a deployed model). Separately, nested CV below re-fits LASSO
# INSIDE each outer fold - not on the full set beforehand - so the outer AUC
# is an honest estimate of the whole pipeline (feature selection + classifier),
# not inflated by selection leakage. Repeated nested CV + a permutation test
# are run for stability and a null baseline. THIS SCRIPT LOCKS ALL THRESHOLDS
# BELOW ONCE RUN - do not re-tune after seeing Stage 6 validation AUCs.
# =============================================================================

library(glmnet)
library(caret)
library(randomForest)
library(xgboost)
library(e1071)
library(pROC)

set.seed(42)

# ---- Config (frozen once this script runs) -----------------------------------
N_REPEATS_MODEL_SELECTION <- 3    # outer repeats while picking RF/SVM/XGB
N_REPEATS_FINAL_ESTIMATE  <- 10   # outer repeats for the reported AUC of the winner
N_PERMUTATIONS            <- 200  # label-shuffle null distribution
LASSO_S                   <- "lambda.1se"

discovery <- readRDS("data/processed/discovery_normalized.rds")
hub_genes <- read.csv("results/tables/hub_genes.csv")$gene
if (length(hub_genes) < 3) {
  hub_genes <- read.csv("results/tables/module_gene_full_table.csv")$gene[1:50]
  cat("Using top 50 module genes instead of strict hub list (too few passed threshold).\n")
}

make_binary_outcome <- function(status) factor(ifelse(status == "ruptured", "ruptured", "not_ruptured"),
                                                levels = c("not_ruptured", "ruptured"))

X_full <- t(discovery$expr[intersect(hub_genes, rownames(discovery$expr)), ])
y_full <- make_binary_outcome(discovery$status)
cat(sprintf("Discovery matrix for ML: %d samples x %d hub genes\n", nrow(X_full), ncol(X_full)))
cat("Class balance:\n"); print(table(y_full))

# ---- LASSO gene selection - REUSABLE, called fresh inside every outer fold --
select_lasso_panel <- function(X_train, y_train) {
  cvfit <- cv.glmnet(X_train, y_train, family = "binomial", alpha = 1, nfolds = 5)
  coefs <- coef(cvfit, s = LASSO_S)
  genes <- setdiff(rownames(coefs)[which(coefs[, 1] != 0)], "(Intercept)")
  if (length(genes) < 2) {
    # lambda.1se degenerated to a near-null model (common with small n,
    # many predictors). Fall back to lambda.min - still fully
    # cross-validated, just less conservative. This is a fixed rule
    # applied identically everywhere this function runs (every fold AND
    # the final fit) - not a post-hoc tune based on results.
    coefs <- coef(cvfit, s = "lambda.min")
    genes <- setdiff(rownames(coefs)[which(coefs[, 1] != 0)], "(Intercept)")
  }
  genes
}

# ---- One outer fold: select on train only, evaluate on held-out test --------
run_outer_fold <- function(train_idx, test_idx, method, X = X_full, y = y_full,
                            fixed_hyperparams = FALSE) {
  train_X <- X[train_idx, , drop = FALSE]; train_y <- y[train_idx]
  test_X  <- X[test_idx, , drop = FALSE];  test_y  <- y[test_idx]

  panel <- select_lasso_panel(train_X, train_y)
  if (length(panel) < 2) return(NA)  # degenerate fold, skip

  train_X_panel <- train_X[, panel, drop = FALSE]
  test_X_panel  <- test_X[, intersect(panel, colnames(test_X)), drop = FALSE]
  missing <- setdiff(panel, colnames(test_X_panel))
  if (length(missing)) {
    pad <- matrix(0, nrow = nrow(test_X_panel), ncol = length(missing),
                   dimnames = list(rownames(test_X_panel), missing))
    test_X_panel <- cbind(test_X_panel, pad)[, panel, drop = FALSE]
  }

  if (fixed_hyperparams) {
    ctrl <- trainControl(method = "none", classProbs = TRUE)
  } else {
    ctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)
  }
  model <- tryCatch(
    train(x = train_X_panel, y = train_y, method = method, trControl = ctrl, metric = "ROC"),
    error = function(e) NULL
  )
  if (is.null(model)) return(NA)

  pred_prob <- predict(model, test_X_panel, type = "prob")[, "ruptured"]
  if (length(unique(test_y)) < 2) return(NA)
  as.numeric(auc(roc(test_y, pred_prob, quiet = TRUE, levels = rev(levels(test_y)))))
}

repeated_nested_cv <- function(method, n_repeats, y = y_full) {
  all_auc <- c()
  for (r in 1:n_repeats) {
    folds <- createFolds(y, k = 5, list = TRUE)
    for (fold in folds) {
      test_idx <- fold
      train_idx <- setdiff(seq_along(y), fold)
      all_auc <- c(all_auc, run_outer_fold(train_idx, test_idx, method))
    }
  }
  all_auc[!is.na(all_auc)]
}

# ---- Step 1: model selection via repeated nested CV, proper leakage-free ----
cat(sprintf("\n>> Model selection: %d-repeat nested CV per method (this takes a while)\n",
            N_REPEATS_MODEL_SELECTION))
cat(">> Random Forest...\n"); rf_auc  <- repeated_nested_cv("rf", N_REPEATS_MODEL_SELECTION)
cat(">> SVM (radial)...\n");  svm_auc <- repeated_nested_cv("svmRadial", N_REPEATS_MODEL_SELECTION)

model_selection_results <- data.frame(
  model = c("RandomForest", "SVM_radial", "XGBoost"),
  mean_AUC = c(mean(rf_auc), mean(svm_auc), mean(xgb_auc)),
  sd_AUC = c(sd(rf_auc), sd(svm_auc), sd(xgb_auc)),
  n_folds = c(length(rf_auc), length(svm_auc), length(xgb_auc))
)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(model_selection_results, "results/tables/model_selection_nested_cv.csv", row.names = FALSE)
cat("\n=== Model selection (nested, leakage-free) ===\n")
print(model_selection_results)

best_method <- model_selection_results$model[which.max(model_selection_results$mean_AUC)]
best_method_caret <- switch(best_method, "RandomForest" = "rf", "SVM_radial" = "svmRadial", "XGBoost" = "xgbTree")
cat(sprintf("\nBest method: %s\n", best_method))

# ---- Step 2: final honest AUC estimate for the winner, more repeats ---------
cat(sprintf("\n>> Final estimate: %d-repeat nested CV for %s\n", N_REPEATS_FINAL_ESTIMATE, best_method))
final_auc_dist <- repeated_nested_cv(best_method_caret, N_REPEATS_FINAL_ESTIMATE)
final_estimate <- data.frame(
  model = best_method,
  mean_AUC = mean(final_auc_dist),
  sd_AUC = sd(final_auc_dist),
  ci_lower = quantile(final_auc_dist, 0.025),
  ci_upper = quantile(final_auc_dist, 0.975),
  n_folds = length(final_auc_dist)
)
write.csv(final_estimate, "results/tables/nested_cv_discovery_auc.csv", row.names = FALSE)
cat("\n=== Final nested CV AUC estimate (report this, not a single-split number) ===\n")
print(final_estimate)

# ---- Step 3: permutation test - null distribution under shuffled labels -----
cat(sprintf("\n>> Permutation test: %d shuffles x single 5-fold nested CV each (fixed hyperparams for speed)\n",
            N_PERMUTATIONS))
perm_aucs <- numeric(N_PERMUTATIONS)
pb_step <- max(1, floor(N_PERMUTATIONS / 20))
for (p in 1:N_PERMUTATIONS) {
  y_shuffled <- sample(y_full)
  folds <- createFolds(y_shuffled, k = 5, list = TRUE)
  fold_aucs <- sapply(folds, function(fold) {
    test_idx <- fold
    train_idx <- setdiff(seq_along(y_shuffled), fold)
    run_outer_fold(train_idx, test_idx, best_method_caret, y = y_shuffled, fixed_hyperparams = TRUE)
  })
  perm_aucs[p] <- mean(fold_aucs, na.rm = TRUE)
  if (p %% pb_step == 0) cat(sprintf("  ...%d/%d permutations done\n", p, N_PERMUTATIONS))
}
perm_aucs <- perm_aucs[!is.na(perm_aucs)]
observed_auc <- final_estimate$mean_AUC
perm_p_value <- mean(perm_aucs >= observed_auc)

write.csv(data.frame(permuted_AUC = perm_aucs), "results/tables/permutation_null_distribution.csv", row.names = FALSE)
cat(sprintf("\n=== Permutation test ===\nObserved AUC: %.3f | Null mean: %.3f (sd %.3f) | p = %.4f\n",
            observed_auc, mean(perm_aucs), sd(perm_aucs), perm_p_value))

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
png("results/figures/permutation_null_distribution.png", width = 600, height = 450)
hist(perm_aucs, breaks = 30, main = "Permutation null distribution vs observed AUC",
     xlab = "AUC", xlim = range(c(perm_aucs, observed_auc)))
abline(v = observed_auc, col = "red", lwd = 2)
legend("topright", legend = sprintf("Observed AUC = %.3f", observed_auc), col = "red", lwd = 2, bty = "n")
dev.off()

# ---- Step 4: train the FINAL deployed model on ALL discovery data -----------
# This is the model/panel that Stage 6/7/8/9 actually use - selected once,
# frozen from here on. Nested CV above is an evaluation of this same pipeline,
# not a different model.
final_panel <- select_lasso_panel(X_full, y_full)
cat(sprintf("\nFinal deployed panel (LASSO on full discovery set): %d genes\n", length(final_panel)))
print(final_panel)
write.csv(data.frame(gene = final_panel), "results/tables/final_gene_panel.csv", row.names = FALSE)

if (length(final_panel) < 2) {
  stop("LASSO selected fewer than 2 genes on the full set - revisit hub gene threshold in 04_wgcna_modules.R.")
}

X_panel_full <- X_full[, final_panel, drop = FALSE]
ctrl_final <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)
final_model <- train(x = X_panel_full, y = y_full, method = best_method_caret, trControl = ctrl_final, metric = "ROC")
dir.create("results/model", showWarnings = FALSE, recursive = TRUE)
saveRDS(final_model, "results/model/final_classifier.rds")

cat("\n=== Stage 5 complete ===\n")
cat(sprintf("Deployed model: %s, %d-gene panel, frozen in results/tables/final_gene_panel.csv\n",
            best_method, length(final_panel)))
cat(sprintf("Honest nested-CV performance estimate: AUC %.3f (95%% CI %.3f-%.3f), permutation p = %.4f\n",
            final_estimate$mean_AUC, final_estimate$ci_lower, final_estimate$ci_upper, perm_p_value))
cat("Hand results/model/final_classifier.rds + final_gene_panel.csv to Role 4 for 06/07.\n")
cat("These thresholds and this panel are now LOCKED - do not re-tune after seeing Stage 6 validation results.\n")

