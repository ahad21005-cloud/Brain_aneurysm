# =============================================================================
# 05b_permutation_test_corrected.R
# Corrected nested-CV permutation test for the 4-gene LASSO classifier.
#
# What this fixes relative to the earlier permutation run:
#   1. EVERY shuffle yields a valid AUC. Earlier, shuffles where the fold-internal
#      LASSO picked < 2 genes were dropped (136/500), so p was conditional on the
#      subset that "worked". Here the gene-selection rule always returns >= 2 genes
#      (lambda.1se -> lambda.min -> smallest-lambda-with->=2-genes on the path ->
#      top-2 univariate), applied IDENTICALLY to real and permuted labels.
#   2. Observed and permuted statistics use the IDENTICAL pipeline (same label-blind
#      5,000-gene variance pool, same LASSO rule, same fixed random forest, same
#      explicit AUC direction). No model selection happens on the real labels, so
#      nothing is tuned on the observed data that the null did not also get.
#   3. p = (k + 1) / (n + 1) over ALL n permutations.
#
# Run from the project root, after 03_normalize_discovery.R has produced
#   data/processed/discovery_normalized.rds   (list: expr = genes x samples, status)
# Needs: glmnet, randomForest, pROC  (parallel is part of base R)
# =============================================================================
suppressPackageStartupMessages({ library(glmnet); library(randomForest); library(pROC); library(parallel) })

# ---- settings -----------------------------------------------------------------
N_PERM          <- 500     # label shuffles
N_REPEATS_NULL  <- 1       # repeats of 5-fold CV per shuffle
N_REPEATS_OBS   <- 20      # repeats of 5-fold CV for the observed statistic
K_FOLDS         <- 5
POOL_N          <- 5000    # label-blind variance pool (same as WGCNA input)
N_TREES         <- 500
SEED            <- 42
N_CORES         <- if (.Platform$OS.type == "windows") 1L else max(1L, detectCores() - 1L)
IN_FILE         <- "data/processed/discovery_normalized.rds"
OUT_DIR         <- "results/tables"
# -------------------------------------------------------------------------------

set.seed(SEED); RNGkind("L'Ecuyer-CMRG"); set.seed(SEED)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

disc   <- readRDS(IN_FILE)
expr   <- disc$expr
y_real <- factor(ifelse(disc$status == "ruptured", "ruptured", "not_ruptured"),
                 levels = c("not_ruptured", "ruptured"))
stopifnot(ncol(expr) == length(y_real))
cat(sprintf("Samples: %d (%d ruptured, %d not)\n", length(y_real), sum(y_real == "ruptured"), sum(y_real != "ruptured")))

# label-blind candidate pool: variance only, computed once, never touches labels
gene_var <- apply(expr, 1, var)
pool     <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(POOL_N, length(gene_var)))]
X        <- t(expr[pool, ])
cat(sprintf("Candidate pool: %d genes (variance-ranked, label-blind)\n", ncol(X)))

nonzero <- function(cv, s) {
  b <- coef(cv, s = s); setdiff(rownames(b)[b[, 1] != 0], "(Intercept)")
}

# gene selection that ALWAYS returns >= 2 genes; also reports which rule fired
select_panel <- function(Xtr, ytr) {
  cv <- cv.glmnet(Xtr, ytr, family = "binomial", alpha = 1, nfolds = 5)
  g <- nonzero(cv, "lambda.1se"); rule <- "lambda.1se"
  if (length(g) < 2) { g <- nonzero(cv, "lambda.min"); rule <- "lambda.min" }
  if (length(g) < 2) {
    k <- which(cv$glmnet.fit$df >= 2)
    if (length(k)) { g <- nonzero(cv, cv$glmnet.fit$lambda[min(k)]); rule <- "path_ge2" }
  }
  if (length(g) < 2) {
    tt <- apply(Xtr, 2, function(x) tryCatch(abs(unname(t.test(x ~ ytr)$statistic)), error = function(e) 0))
    g <- names(sort(tt, decreasing = TRUE))[1:2]; rule <- "top2_univariate"
  }
  list(genes = g, rule = rule)
}

# stratified fold assignment
make_folds <- function(y, k) {
  f <- integer(length(y))
  for (cl in levels(y)) { i <- which(y == cl); f[i] <- sample(rep_len(seq_len(k), length(i))) }
  f
}

fold_auc <- function(tr, te, y) {
  sel <- select_panel(X[tr, , drop = FALSE], y[tr])
  rf  <- randomForest(x = X[tr, sel$genes, drop = FALSE], y = y[tr], ntree = N_TREES)
  p   <- predict(rf, X[te, sel$genes, drop = FALSE], type = "prob")[, "ruptured"]
  a   <- as.numeric(auc(roc(y[te], p, levels = c("not_ruptured", "ruptured"), direction = "<", quiet = TRUE)))
  list(auc = a, rule = sel$rule)
}

cv_stat <- function(y, n_rep) {
  aucs <- c(); rules <- c()
  for (r in seq_len(n_rep)) {
    f <- make_folds(y, K_FOLDS)
    for (k in seq_len(K_FOLDS)) {
      te <- which(f == k); tr <- which(f != k)
      o <- fold_auc(tr, te, y); aucs <- c(aucs, o$auc); rules <- c(rules, o$rule)
    }
  }
  list(auc = mean(aucs), fold_aucs = aucs, rules = rules)
}

# ---- observed -------------------------------------------------------------------
cat(sprintf("\nObserved statistic: %d repeats x %d folds ...\n", N_REPEATS_OBS, K_FOLDS))
t0  <- Sys.time()
obs <- cv_stat(y_real, N_REPEATS_OBS)
cat(sprintf("Observed mean fold AUC = %.3f (SD of folds %.3f)  [%.1f min]\n",
            obs$auc, sd(obs$fold_aucs), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("Selection rule used in observed folds:\n"); print(round(prop.table(table(obs$rules)), 3))

# ---- permutations ---------------------------------------------------------------
cat(sprintf("\nPermutations: %d shuffles x %d repeat(s) x %d folds, %d core(s) ...\n",
            N_PERM, N_REPEATS_NULL, K_FOLDS, N_CORES))
one_perm <- function(i) {
  o <- cv_stat(sample(y_real), N_REPEATS_NULL)
  list(auc = o$auc, rules = o$rules)
}
res <- mclapply(seq_len(N_PERM), one_perm, mc.cores = N_CORES, mc.set.seed = TRUE)
ok  <- vapply(res, function(r) is.list(r) && is.finite(r$auc), logical(1))
cat(sprintf("Valid permutations: %d / %d\n", sum(ok), N_PERM))
null_auc  <- vapply(res[ok], function(r) r$auc, numeric(1))
null_rule <- unlist(lapply(res[ok], function(r) r$rules))

k <- sum(null_auc >= obs$auc)
p <- (k + 1) / (length(null_auc) + 1)
cat(sprintf("\n=== RESULT ===\nObserved AUC %.3f | null mean %.3f (SD %.3f) | exceedances %d/%d | p = (k+1)/(n+1) = %.4f\n",
            obs$auc, mean(null_auc), sd(null_auc), k, length(null_auc), p))
cat("Selection rule used in permuted folds:\n"); print(round(prop.table(table(null_rule)), 3))

write.csv(data.frame(permuted_AUC = null_auc), file.path(OUT_DIR, "permutation_null_corrected_v2.csv"), row.names = FALSE)
write.csv(data.frame(observed_AUC = obs$auc, sd_folds = sd(obs$fold_aucs),
                     pct_lo = quantile(obs$fold_aucs, .025), pct_hi = quantile(obs$fold_aucs, .975),
                     null_mean = mean(null_auc), null_sd = sd(null_auc),
                     n_perm_valid = length(null_auc), exceed = k, p_value = p),
          file.path(OUT_DIR, "permutation_summary_corrected_v2.csv"), row.names = FALSE)
png("results/figures/permutation_null_corrected_v2.png", width = 1800, height = 1350, res = 300)
hist(null_auc, breaks = 30, col = "grey85", border = "grey40", main = "Permutation null vs observed AUC",
     xlab = "AUC", xlim = range(c(null_auc, obs$auc)))
abline(v = obs$auc, col = "red", lwd = 2)
legend("topleft", legend = sprintf("Observed AUC = %.3f (p = %.3f)", obs$auc, p), col = "red", lwd = 2, bty = "n")
dev.off()
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo_permutation.txt"))
cat("\nDone. Send me the === RESULT === block and the two 'Selection rule' tables.\n")
