# =============================================================================
# 09_scRNA_deconvolution.R  (Role 5 owns this)
#
# SCOPE NOTE (state this in the writeup): GSE256490 contains no IA samples
# and no cell-type labels. This script uses a small subset of the LEAST
# diseased samples (unsorted fetal CNS + epilepsy-resection "normal-adjacent"
# cortex, not tumor/AVM/metastasis samples) to build a COARSE, GENERAL
# vascular/immune reference (endothelial, SMC/pericyte, immune, fibroblast,
# other) via standard clustering + canonical marker-gene scoring - not
# full manual expert annotation, and not IA-specific. This reference is only
# valid for describing which broad cell type the panel genes are enriched
# in - not for any IA-specific mechanistic claim.
# =============================================================================

library(Seurat)
library(Matrix)
library(nnls)

raw_dir <- "data/raw/GSE256490/extracted"

# ---- Samples chosen: unsorted, non-diseased-adjacent tissue only ------------
sample_prefixes <- c(
  "GSM8101564_Fetal_CNS_20-10-2020_unsorted_3pr_v3",
  "GSM8101569_TL_10-10-2019_unsorted",
  "GSM8101572_TLE_White-matter_Unsorted_3pr_v3",
  "GSM8101573_TLE_Grey-matter_Unsorted_3pr_v3"
)

load_one_sample <- function(prefix, dir) {
  files <- list.files(dir, pattern = paste0("^", prefix), full.names = TRUE)
  barcode_f <- grep("barcodes", files, value = TRUE)
  feature_f <- grep("features|genes", files, value = TRUE)
  matrix_f  <- grep("matrix", files, value = TRUE)
  stopifnot(length(barcode_f) == 1, length(feature_f) == 1, length(matrix_f) == 1)

  mat <- readMM(matrix_f)
  barcodes <- readLines(barcode_f)
  features <- read.delim(feature_f, header = FALSE)
  rownames(mat) <- make.unique(features[[2]])  # gene symbols, col 2 of features.tsv
  colnames(mat) <- barcodes
  colnames(mat) <- paste0(prefix, "_", colnames(mat))  # avoid barcode collisions across samples
  mat
}

cat(">> Loading samples...\n")
mats <- lapply(sample_prefixes, load_one_sample, dir = raw_dir)
common_genes <- Reduce(intersect, lapply(mats, rownames))
mats <- lapply(mats, function(m) m[common_genes, ])
combined <- do.call(cbind, mats)
cat(sprintf("Combined: %d genes x %d cells across %d samples\n",
            nrow(combined), ncol(combined), length(sample_prefixes)))

# ---- Seurat pipeline: QC -> normalize -> cluster -----------------------------
so <- CreateSeuratObject(counts = combined, min.cells = 3, min.features = 200)
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so <- subset(so, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
cat(sprintf("After QC filtering: %d cells remain\n", ncol(so)))

so <- NormalizeData(so)
so <- FindVariableFeatures(so, nfeatures = 2000)
so <- ScaleData(so, features = VariableFeatures(so))
so <- RunPCA(so, npcs = 30, verbose = FALSE)
so <- FindNeighbors(so, dims = 1:30)
so <- FindClusters(so, resolution = 0.5)

cat(sprintf("Clusters found: %d\n", length(unique(Idents(so)))))

# ---- Marker-based cluster annotation -----------------------------------------
marker_sets <- list(
  Endothelial = c("PECAM1", "CDH5", "VWF", "CLDN5"),
  SMC_Pericyte = c("ACTA2", "MYH11", "PDGFRB", "RGS5", "TAGLN"),
  Immune = c("PTPRC", "CD68", "CD3E", "AIF1"),
  Fibroblast = c("COL1A1", "DCN", "LUM", "PDGFRA"),
  Neuronal_Glial = c("RBFOX3", "GFAP", "OLIG1", "SNAP25")
)

for (ct in names(marker_sets)) {
  genes_present <- intersect(marker_sets[[ct]], rownames(so))
  if (length(genes_present) == 0) next
  so <- AddModuleScore(so, features = list(genes_present), name = ct, seed = 42)
}
score_cols <- grep(paste0("^(", paste(names(marker_sets), collapse = "|"), ")\\d*$"),
                    colnames(so@meta.data), value = TRUE)

cluster_scores <- aggregate(so@meta.data[, score_cols], by = list(cluster = Idents(so)), mean)
cluster_assignment <- names(marker_sets)[apply(cluster_scores[, score_cols], 1, which.max)]
names(cluster_assignment) <- cluster_scores$cluster

cat("\nCluster -> cell type assignment:\n")
print(data.frame(cluster = cluster_scores$cluster, assigned = cluster_assignment))

so$cell_type <- unname(cluster_assignment[as.character(Idents(so))])

# ---- Build signature matrix: mean expression per assigned cell type ---------
expr_data <- GetAssayData(so, layer = "data")
sig_matrix <- sapply(split(colnames(so), so$cell_type),
                      function(cells) rowMeans(as.matrix(expr_data[, cells, drop = FALSE])))
cat(sprintf("\nSignature matrix: %d genes x %d cell types\n", nrow(sig_matrix), ncol(sig_matrix)))
print(colnames(sig_matrix))

dir.create("results/model", showWarnings = FALSE, recursive = TRUE)
saveRDS(sig_matrix, "results/model/sc_signature_matrix.rds")
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(data.frame(cluster = cluster_scores$cluster, assigned = cluster_assignment),
          "results/tables/scRNA_cluster_annotation.csv", row.names = FALSE)

# =============================================================================
# STEP 2: NNLS deconvolution of each bulk cohort against sig_matrix
# =============================================================================

deconvolve_cohort <- function(bulk_expr, sig_matrix) {
  common <- intersect(rownames(bulk_expr), rownames(sig_matrix))
  cat(sprintf("  %d genes shared between bulk cohort and signature matrix\n", length(common)))
  B <- as.matrix(bulk_expr[common, , drop = FALSE])
  S <- as.matrix(sig_matrix[common, , drop = FALSE])
  props <- t(apply(B, 2, function(col) { fit <- nnls(S, col); p <- fit$x; p / sum(p) }))
  colnames(props) <- colnames(S)
  rownames(props) <- colnames(B)
  as.data.frame(props)
}

cohorts_to_deconvolve <- list(
  discovery = "data/processed/discovery_normalized.rds",
  GSE13353  = "data/processed/GSE13353_normalized.rds",
  GSE15629  = "data/processed/GSE15629_normalized.rds",
  GSE54083  = "data/processed/GSE54083_normalized.rds"
)

all_proportions <- list()
for (nm in names(cohorts_to_deconvolve)) {
  cat(sprintf("\n>> Deconvolving %s\n", nm))
  d <- readRDS(cohorts_to_deconvolve[[nm]])
  props <- deconvolve_cohort(d$expr, sig_matrix)
  props$status <- d$status
  all_proportions[[nm]] <- props
  write.csv(props, sprintf("results/tables/celltype_proportions_%s.csv", nm), row.names = TRUE)
}

# =============================================================================
# STEP 3: correlate panel genes against cell-type proportions, FDR-corrected
# =============================================================================

selected_genes <- read.csv("results/tables/final_gene_panel.csv")$gene
discovery <- readRDS("data/processed/discovery_normalized.rds")
disc_props <- all_proportions[["discovery"]]

cor_results <- do.call(rbind, lapply(selected_genes, function(g) {
  if (!(g %in% rownames(discovery$expr))) return(NULL)
  gene_expr <- discovery$expr[g, ]
  do.call(rbind, lapply(setdiff(colnames(disc_props), "status"), function(ct) {
    test <- cor.test(gene_expr, disc_props[[ct]], method = "spearman")
    data.frame(gene = g, cell_type = ct, rho = unname(test$estimate), p_value = test$p.value)
  }))
}))
cor_results$p_adj_FDR <- p.adjust(cor_results$p_value, method = "BH")
cor_results <- cor_results[order(cor_results$p_adj_FDR), ]
write.csv(cor_results, "results/tables/gene_celltype_correlations.csv", row.names = FALSE)

cat("\n=== Stage 8 complete ===\n")
cat(sprintf("Significant gene-celltype associations (FDR < 0.05): %d / %d\n",
            sum(cor_results$p_adj_FDR < 0.05, na.rm = TRUE), nrow(cor_results)))
cat("\nREMEMBER for writeup: reference built from 4 non-diseased-adjacent samples\n")
cat("of a general CNS atlas (no IA samples, no ground-truth cell labels) using\n")
cat("marker-score cluster annotation - report as a coarse/general reference,\n")
cat("not an IA-specific one.\n")

