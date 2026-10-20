# =============================================================================
# 06_cross_platform_transfer.R - apply the FROZEN model, unchanged, to the
# independent cohorts and run the expression-range diagnostic.
#   Tissue (ruptured vs unruptured): GSE13353, GSE15629, GSE54083
#   Blood (SAH vs control, exploratory): GSE36791
# Reports AUC (direction "<", bootstrap percentile CI), the per-gene raw ranges
# vs discovery (e.g. LIMCH1 4.06-10.60 in discovery vs 15.22-16.52 in GSE54083)
# and the spread of predicted probabilities (near-constant => scale mismatch).
# Nothing is retrained or recalibrated here.
# =============================================================================
suppressPackageStartupMessages({ library(randomForest); library(pROC) })
set.seed(42)
fm <- readRDS("results/model/final_classifier.rds"); genes <- fm$genes
cohorts <- list(GSE13353 = c("ruptured", "unruptured"), GSE15629 = c("ruptured", "unruptured"),
                GSE54083 = c("ruptured", "unruptured"), GSE36791 = c("SAH", "control"))
disc <- readRDS("data/processed/discovery_normalized.rds")

boot_auc <- function(y, p, B = 2000) {
  a <- replicate(B, { i <- sample(seq_along(y), replace = TRUE)
    if (length(unique(y[i])) < 2) NA else as.numeric(auc(roc(y[i], p[i], direction = "<", quiet = TRUE))) })
  quantile(a, c(.025, .975), na.rm = TRUE)
}
res <- list(); rng <- list()
for (nm in names(cohorts)) {
  d <- readRDS(sprintf("data/processed/%s_normalized.rds", nm))
  pos <- cohorts[[nm]][1]; neg <- cohorts[[nm]][2]
  keep <- d$status %in% c(pos, neg); st <- d$status[keep]; E <- d$expr[, keep, drop = FALSE]
  present <- intersect(genes, rownames(E))
  for (g in genes) rng[[length(rng) + 1]] <- data.frame(
    cohort = nm, gene = g, in_cohort = g %in% rownames(E),
    disc_min = fm$disc_range[g, 1], disc_max = fm$disc_range[g, 2],
    val_min = if (g %in% rownames(E)) min(E[g, ]) else NA, val_max = if (g %in% rownames(E)) max(E[g, ]) else NA)
  if (length(present) < length(genes)) { cat(sprintf("[%s] missing panel genes: %s - skipped\n", nm, paste(setdiff(genes, present), collapse = ","))); next }
  X <- t(E[genes, , drop = FALSE]); colnames(X) <- genes
  p <- predict(fm$model, X, type = "prob")[, "ruptured"]
  y <- factor(ifelse(st == pos, 1, 0), levels = c(0, 1))
  a <- as.numeric(auc(roc(y, p, direction = "<", quiet = TRUE))); ci <- boot_auc(y, p)
  res[[nm]] <- data.frame(cohort = nm, n_pos = sum(y == 1), n_neg = sum(y == 0), AUC = a,
                          ci_lo = ci[1], ci_hi = ci[2], pred_sd = sd(p), pred_unique = length(unique(round(p, 4))))
}
res <- do.call(rbind, res); rng <- do.call(rbind, rng); rownames(res) <- NULL; rownames(rng) <- NULL
write.csv(res, "results/tables/frozen_model_transfer.csv", row.names = FALSE)
write.csv(rng, "results/tables/cross_platform_range_diagnostic.csv", row.names = FALSE)
cat("\n=== Frozen-model transfer ===\n"); print(res, digits = 3)
cat("\n=== Raw expression ranges, discovery vs cohort ===\n"); print(rng, digits = 3)
cat("\nNon-overlapping ranges + near-constant predictions = scale mismatch, not biology.\n")
