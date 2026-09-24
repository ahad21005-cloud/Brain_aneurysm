# =============================================================================
# 03_normalize_discovery.R  (Role 1 owns this)
# Normalizes GSE122897 ONLY. No batch correction needed - this is a single
# RNA-seq batch carrying all three classes (ruptured/unruptured/control),
# so there is no second cohort to merge against.
#
# NOTE: GSE122897 series matrix ships with NO expression data (0 genes) -
# counts live in the supplementary file GSE122897_readCounts_raw.txt.gz
# instead. pData(eset) is matched to the counts file by sample title
# (CO1/IA1-style), NOT by row order.
# =============================================================================

library(Biobase)
library(DESeq2)
source("scripts/utils/label_utils.R")
source("scripts/utils/gene_mapping.R")

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

eset <- readRDS("data/raw/GSE122897_eset.rds")

counts <- read.delim("data/raw/GSE122897/GSE122897_readCounts_raw.txt.gz", row.names = 1)
counts <- as.matrix(counts)

cat("Feature ID format check (first 10 row names):\n")
print(head(rownames(counts), 10))

pd <- pData(eset)
stopifnot(all(colnames(counts) %in% pd$title))
pd_aligned <- pd[match(colnames(counts), pd$title), ]
stopifnot(identical(pd_aligned$title, colnames(counts)))

symbol_map <- map_rnaseq_features(rownames(counts))

counts_int <- round(counts)
mode(counts_int) <- "integer"

status <- label_gse122897(pd_aligned)

n_dropped <- sum(is.na(status))
if (n_dropped > 0) {
  cat(sprintf("Dropping %d sample(s) with unresolved rupture status.\n", n_dropped))
  counts_int <- counts_int[, !is.na(status), drop = FALSE]
  status <- status[!is.na(status)]
}
check_labels(status, "GSE122897")

col_data <- data.frame(status = factor(status, levels = c("control", "unruptured", "ruptured")))
rownames(col_data) <- colnames(counts_int)

dds <- DESeqDataSetFromMatrix(countData = counts_int, colData = col_data, design = ~status)
dds <- dds[rowSums(counts(dds)) > 10, ]
vsd <- vst(dds, blind = FALSE)
expr_vst <- assay(vsd)

expr_gene <- collapse_to_gene_symbol(expr_vst, symbol_map)

out <- list(expr = expr_gene, status = status, sample_id = colnames(expr_gene))
saveRDS(out, "data/processed/discovery_normalized.rds")

cat(sprintf("\nSaved: data/processed/discovery_normalized.rds\n"))
cat(sprintf("Discovery matrix: %d genes x %d samples (%d dropped)\n",
            nrow(expr_gene), ncol(expr_gene), n_dropped))
cat("=== Stage 3 (discovery) complete. ===\n")

