# =============================================================================
# 04_wgcna_modules.R  (Role 3 owns this)
# Runs directly on the discovery cohort (GSE122897) - no merge, no batch
# correction needed, since it is a single RNA-seq batch with all three
# classes already. Can also submit via hpc/slurm_wgcna.sh if too slow locally.
# =============================================================================

library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

dat <- readRDS("data/processed/discovery_normalized.rds")
expr <- t(dat$expr)
status <- dat$status

gene_var <- apply(expr, 2, var, na.rm = TRUE)
keep_n <- min(5000, ncol(expr))
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:keep_n]
expr_filt <- expr[, top_genes]

gsg <- goodSamplesGenes(expr_filt, verbose = 3)
if (!gsg$allOK) {
  expr_filt <- expr_filt[gsg$goodSamples, gsg$goodGenes]
  status <- status[gsg$goodSamples]
}

powers <- 1:20
sft <- pickSoftThreshold(expr_filt, powerVector = powers, verbose = 2)

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
png("results/figures/wgcna_soft_threshold.png", width = 900, height = 450)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology R^2", main = "Scale independence")
abline(h = 0.85, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", main = "Mean connectivity")
dev.off()

chosen_power <- sft$powerEstimate
if (is.na(chosen_power)) {
  cat("WARNING: no power reached R^2 > 0.85. Inspect the saved plot and set manually.\n")
  chosen_power <- 6
}
cat(sprintf("Using soft-thresholding power: %d\n", chosen_power))

net <- blockwiseModules(
  expr_filt, power = chosen_power, networkType = "signed", TOMType = "signed",
  minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, saveTOMs = TRUE,
  saveTOMFileBase = "results/model/discovery_TOM", verbose = 3, maxBlockSize = 6000
)

status_numeric <- as.numeric(factor(status, levels = c("control", "unruptured", "ruptured")))
MEs <- net$MEs
module_trait_cor <- cor(MEs, status_numeric, use = "p")
module_trait_p <- corPvalueStudent(module_trait_cor, nSamples = nrow(expr_filt))

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
module_trait_df <- data.frame(module = colnames(MEs),
                               correlation = module_trait_cor[, 1],
                               p_value = module_trait_p[, 1])
module_trait_df <- module_trait_df[order(module_trait_df$p_value), ]
write.csv(module_trait_df, "results/tables/module_trait_relationships.csv", row.names = FALSE)

cat("\n=== Module-trait relationships (sorted by p-value) ===\n")
print(head(module_trait_df, 10))

top_module_label <- as.numeric(gsub("ME", "", module_trait_df$module[1]))
gene_module_membership <- as.data.frame(cor(expr_filt, MEs, use = "p"))
gene_significance <- as.numeric(cor(expr_filt, status_numeric, use = "p"))

top_module_genes <- colnames(expr_filt)[net$colors == top_module_label]
mm_col <- paste0("ME", top_module_label)

hub_table <- data.frame(
  gene = top_module_genes,
  module_membership = gene_module_membership[top_module_genes, mm_col],
  gene_significance = gene_significance[match(top_module_genes, colnames(expr_filt))]
)
hub_table <- hub_table[order(-abs(hub_table$module_membership)), ]
hub_genes <- hub_table[abs(hub_table$module_membership) > 0.8 & abs(hub_table$gene_significance) > 0.2, ]

write.csv(hub_table, "results/tables/module_gene_full_table.csv", row.names = FALSE)
write.csv(hub_genes, "results/tables/hub_genes.csv", row.names = FALSE)

cat(sprintf("\nTop rupture-associated module: %s (n=%d genes)\n", module_trait_df$module[1], length(top_module_genes)))
cat(sprintf("Hub genes extracted (|MM|>0.8, |GS|>0.2): %d\n", nrow(hub_genes)))

saveRDS(net, "results/model/wgcna_net.rds")

cat("\n=== Stage 4 complete. Hand hub_genes.csv to Stage 5 (ML classifier). ===\n")
if (nrow(hub_genes) < 5) {
  cat("WARNING: fewer than 5 hub genes passed threshold - consider loosening\n")
  cat("|MM|/|GS| cutoffs, or use module_gene_full_table.csv directly as LASSO input.\n")
}

