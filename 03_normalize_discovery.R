# =============================================================================
# 03_normalize_discovery.R - GSE122897 (RNA-seq discovery). Single batch, no
# batch correction, never merged with another cohort.
# Counts come from the SUPPLEMENTARY file (series matrix has 0 genes) and are
# matched to GEO metadata by sample TITLE, not by row order.
# Filter: total count > 10; transform: DESeq2 VST (blind = FALSE).
# Writes a SAMPLE AUDIT table - use it to explain 60 GEO samples -> 58 analysed.
# =============================================================================
suppressPackageStartupMessages({ library(Biobase); library(DESeq2) })
source("scripts/utils/label_utils.R"); source("scripts/utils/gene_mapping.R")
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

eset   <- readRDS("data/raw/GSE122897_eset.rds")
cf     <- list.files("data/raw/GSE122897", pattern = "readCounts.*txt(\\.gz)?$", full.names = TRUE)
stopifnot(length(cf) == 1)
counts <- as.matrix(read.delim(cf, row.names = 1, check.names = FALSE))
cat("First feature IDs:\n"); print(head(rownames(counts)))

pd <- pData(eset)
cat(sprintf("GEO samples: %d | count-matrix columns: %d\n", nrow(pd), ncol(counts)))
in_both <- intersect(colnames(counts), pd$title)
stopifnot(length(in_both) > 0)
pd_al  <- pd[match(colnames(counts), pd$title), ]
status <- label_gse122897(pd_al)

# ---- audit: every GEO sample, its label, and why it was (not) analysed ------
audit <- data.frame(geo_accession = rownames(pd), title = pd$title,
                    in_count_matrix = pd$title %in% colnames(counts), stringsAsFactors = FALSE)
audit$label  <- status[match(audit$title, pd_al$title)]
audit$reason <- ifelse(!audit$in_count_matrix, "no column in count matrix",
                ifelse(is.na(audit$label), "rupture status unresolved", "analysed"))
write.csv(audit, "results/tables/GSE122897_sample_audit.csv", row.names = FALSE)
cat("\nSample audit (this explains GEO n vs analysed n):\n"); print(table(audit$reason))

keep <- !is.na(status)
counts_int <- round(counts)[, keep, drop = FALSE]; mode(counts_int) <- "integer"
status <- status[keep]
check_labels(status, "GSE122897")

coldata <- data.frame(status = factor(status, levels = c("control", "unruptured", "ruptured")),
                      row.names = colnames(counts_int))
dds <- DESeqDataSetFromMatrix(counts_int, coldata, design = ~status)
dds <- dds[rowSums(counts(dds)) > 10, ]
vsd <- vst(dds, blind = FALSE)
expr_gene <- collapse_to_gene_symbol(assay(vsd), map_rnaseq_features(rownames(vsd)))

saveRDS(list(expr = expr_gene, status = status, sample_id = colnames(expr_gene)),
        "data/processed/discovery_normalized.rds")
cat(sprintf("\nSaved discovery matrix: %d genes x %d samples\n", nrow(expr_gene), ncol(expr_gene)))
cat("Paper expects 58 samples (16 control / 21 unruptured / 21 ruptured).\n")
