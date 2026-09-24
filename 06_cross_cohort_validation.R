# =============================================================================
# 06_cross_cohort_validation.R  (Role 4 owns this)
# Scores the FROZEN classifier (from 05_ml_classifier.R) independently on
# every validation cohort. Never retrains, never merges cohorts together -
# each one is scored on its own and gets its own AUC + bootstrapped CI.
# =============================================================================

library(caret)
library(pROC)

final_model <- readRDS("results/model/final_classifier.rds")
selected_genes <- read.csv("results/tables/final_gene_panel.csv")$gene

make_binary_outcome <- function(status) factor(ifelse(status == "ruptured", "ruptured", "not_ruptured"),
                                               levels = c("not_ruptured", "ruptured"))

bootstrap_auc_ci <- function(true_labels, pred_probs, n_boot = 2000) {
  aucs <- replicate(n_boot, {
    idx <- sample(seq_along(true_labels), replace = TRUE)
    if (length(unique(true_labels[idx])) < 2) return(NA)
    as.numeric(auc(roc(true_labels[idx], pred_probs[idx], quiet = TRUE, levels = rev(levels(true_labels)))))
  })
  aucs <- aucs[!is.na(aucs)]
  unname(quantile(aucs, c(0.025, 0.5, 0.975)))  # unname() strips the "2.5%"/"50%"/"97.5%" labels at the source
}

validate_cohort <- function(rds_path, cohort_name) {
  d <- readRDS(rds_path)
  
  if (!("ruptured" %in% d$status)) {
    cat(sprintf("[%s] SKIPPED for rupture-AUC purposes - no ruptured samples in this cohort\n", cohort_name))
    return(NULL)
  }
  
  common <- intersect(selected_genes, rownames(d$expr))
  cat(sprintf("\n[%s] %d/%d panel genes present\n", cohort_name, length(common), length(selected_genes)))
  if (length(common) < 2) {
    cat(sprintf("[%s] SKIPPED - too few panel genes present to score reliably.\n", cohort_name))
    return(NULL)
  }
  
  X_val <- t(d$expr[common, , drop = FALSE])
  missing <- setdiff(selected_genes, common)
  if (length(missing)) {
    pad <- matrix(0, nrow = nrow(X_val), ncol = length(missing), dimnames = list(rownames(X_val), missing))
    X_val <- cbind(X_val, pad)[, selected_genes, drop = FALSE]
  }
  
  y_val <- make_binary_outcome(d$status)
  pred_prob <- predict(final_model, X_val, type = "prob")[, "ruptured"]
  ci <- bootstrap_auc_ci(y_val, pred_prob)
  cat(sprintf("[%s] AUC = %.3f (95%% CI: %.3f - %.3f)\n", cohort_name, ci[2], ci[1], ci[3]))
  data.frame(cohort = cohort_name, n_samples = nrow(X_val), AUC = ci[2],
             CI_lower = ci[1], CI_upper = ci[3], row.names = NULL)
}

safe_validate <- function(path, name) {
  tryCatch(validate_cohort(path, name), error = function(e) {
    cat(sprintf("[%s] SKIPPED - file not available yet (%s)\n", name, conditionMessage(e)))
    NULL
  })
}

results <- rbind(
  safe_validate("data/processed/GSE13353_normalized.rds", "GSE13353"),
  safe_validate("data/processed/GSE15629_normalized.rds", "GSE15629"),
  safe_validate("data/processed/GSE54083_normalized.rds", "GSE54083"),
  safe_validate("data/processed/GSE6551_normalized.rds",  "GSE6551")
)

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(results, "results/tables/cross_cohort_validation_auc.csv", row.names = FALSE)

cat("\n=== Cross-cohort validation summary ===\n")
print(results)
cat("\nReport these AUCs WITH their confidence intervals.\n")
cat("\n=== Stage 6 complete. Proceed to 07_blood_validation.R. ===\n")