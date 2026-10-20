# =============================================================================
# 05_ml_classifier.R - CORRECTED classifier (replaces the earlier leaky version).
#   * candidate genes = label-blind top-5,000 variance pool (NOT the WGCNA hub
#     genes: hubs were chosen using the labels, so using them inside CV leaks).
#   * LASSO is re-fitted inside every training fold (select_panel in
#     utils/classifier_utils.R); random forest; AUC with explicit direction.
#   * Observed statistic = mean fold AUC over 20 x 5-fold repeats.
#   * Then the panel is FROZEN: LASSO on all discovery samples -> RF.
# Paper: AUC 0.728 (SD 0.158), permutation p = 0.020 (run 05b), panel
# LIMCH1, ANKRD29, PLCL1, SHF. Exact numbers can differ slightly with R /
# package versions - report what YOUR run gives and keep the seed.
# =============================================================================
source("scripts/utils/classifier_utils.R")
RNGkind("L'Ecuyer-CMRG"); set.seed(42)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/model",  showWarnings = FALSE, recursive = TRUE)

disc <- readRDS("data/processed/discovery_normalized.rds")
y    <- make_outcome(disc$status)
pool <- variance_pool(disc$expr)
X    <- t(disc$expr[pool, ])
cat(sprintf("%d samples (%d ruptured / %d not); candidate pool %d genes\n",
            length(y), sum(y == "ruptured"), sum(y != "ruptured"), ncol(X)))

N_REPEATS <- as.integer(Sys.getenv("IA_REPEATS", "20"))   # 20 in the paper; env var only for quick smoke tests
obs <- cv_stat(X, y, N_REPEATS)
cat(sprintf("Nested-CV mean fold AUC = %.3f (SD %.3f); percentile interval %.3f-%.3f\n",
            obs$auc, sd(obs$fold_aucs), quantile(obs$fold_aucs, .025), quantile(obs$fold_aucs, .975)))
print(round(prop.table(table(obs$rules)), 3))
write.csv(data.frame(mean_AUC = obs$auc, sd_folds = sd(obs$fold_aucs),
                     pct_lo = quantile(obs$fold_aucs, .025), pct_hi = quantile(obs$fold_aucs, .975),
                     n_fold_estimates = length(obs$fold_aucs)),
          "results/tables/nested_cv_discovery_auc.csv", row.names = FALSE)

# ---- freeze the deployed panel + model (all discovery data) -----------------
sel <- select_panel(X, y)
cat("\nFROZEN PANEL:", paste(sel$genes, collapse = ", "), " (rule:", sel$rule, ")\n")
paper_panel <- c("LIMCH1", "ANKRD29", "PLCL1", "SHF")
if (!setequal(sel$genes, paper_panel))
  warning("Panel differs from the paper's (", paste(paper_panel, collapse = ", "),
          "). Check package versions / seed before changing the manuscript.")
rf <- randomForest(x = X[, sel$genes, drop = FALSE], y = y, ntree = N_TREES)
write.csv(data.frame(gene = sel$genes), "results/tables/final_gene_panel.csv", row.names = FALSE)
saveRDS(list(model = rf, genes = sel$genes, pool = pool,
             disc_range = t(apply(disc$expr[sel$genes, , drop = FALSE], 1, range))),
        "results/model/final_classifier.rds")
cat("Saved results/model/final_classifier.rds. THE PANEL IS NOW FROZEN - do not re-tune after seeing validation results.\n")
writeLines(capture.output(sessionInfo()), "results/tables/sessionInfo_classifier.txt")
