# =============================================================================
# 05c_stability_selection.R - bootstrap stability selection of the LASSO panel.
# 100 stratified bootstrap resamples of the discovery samples; LASSO selection
# (same rule as everywhere) on the label-blind pool; report how often each gene
# is selected. Paper: no gene reached >= 60%; best ANKRD29 38%.
# =============================================================================
source("scripts/utils/classifier_utils.R")
RNGkind("L'Ecuyer-CMRG"); set.seed(42)
B <- as.integer(Sys.getenv("IA_BOOT", "100"))   # 100 in the paper
disc <- readRDS("data/processed/discovery_normalized.rds")
y <- make_outcome(disc$status); X <- t(disc$expr[variance_pool(disc$expr), ])
sel_list <- vector("list", B)
for (b in seq_len(B)) {
  idx <- unlist(lapply(levels(y), function(cl) sample(which(y == cl), replace = TRUE)))
  sel_list[[b]] <- select_panel(X[idx, , drop = FALSE], y[idx])$genes
  if (b %% 10 == 0) cat(b, "/", B, "\n")
}
freq <- sort(table(unlist(sel_list)) / B, decreasing = TRUE)
out  <- data.frame(gene = names(freq), selection_frequency = as.numeric(freq))
write.csv(out, "results/tables/stability_selection.csv", row.names = FALSE)
cat("\nTop genes by selection frequency:\n"); print(head(out, 15))
panel <- read.csv("results/tables/final_gene_panel.csv")$gene
cat("\nFrozen panel genes:\n"); print(out[out$gene %in% panel, ])
cat(sprintf("Genes selected in >= 60%% of resamples: %d\n", sum(out$selection_frequency >= 0.6)))
