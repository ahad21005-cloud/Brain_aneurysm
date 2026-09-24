# =============================================================================
# 00_setup.R
# Run ONCE. Installs every R package used across the redesigned pipeline.
# =============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c(
  "GEOquery",                        # dataset download
  "DESeq2",                          # RNA-seq normalization (GSE122897, discovery)
  "limma",                           # microarray normalization
  "sva",                             # kept for reference; NOT used for cross-cohort merging
                                      # in this design (no batch correction step remains -
                                      # discovery is a single cohort), but useful if you
                                      # ever want to sanity-check platform effects
  "WGCNA",                           # co-expression network analysis
  "BisqueRNA",                       # bulk deconvolution using scRNA reference
  "org.Hs.eg.db",                    # Ensembl/Entrez -> gene symbol mapping (GSE122897 RNA-seq)
  "hgu133plus2.db",                  # GPL570 (GSE75436, GSE13353, GSE26969, GSE6551)
  "hugene10sttranscriptcluster.db",  # GPL6244 (GSE15629)
  "illuminaHumanv4.db"                # GPL10558 (GSE36791, blood)
)
BiocManager::install(bioc_pkgs, update = FALSE, ask = FALSE)

cran_pkgs <- c("glmnet", "randomForest", "xgboost", "e1071", "caret", "pROC",
               "tidyverse", "readr", "matrixStats", "Seurat")

new_cran <- cran_pkgs[!(cran_pkgs %in% installed.packages()[, "Package"])]
if (length(new_cran)) install.packages(new_cran)

cat("\n=== Setup check ===\n")
for (pkg in c(bioc_pkgs, cran_pkgs)) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  cat(sprintf("%-35s %s\n", pkg, ifelse(ok, "OK", "MISSING - reinstall manually")))
}

cat("\nNOTES:\n")
cat("- GSE54083 is Agilent (GPL4133): annotation is usually already embedded in\n")
cat("  fData() of the series matrix - 03_normalize_validation.R checks this directly.\n")
cat("- GSE122897 is RNA-seq: check whether its supplementary count matrix uses\n")
cat("  gene symbols or Ensembl IDs before running 03_normalize_discovery.R -\n")
cat("  org.Hs.eg.db handles the conversion if needed.\n")
cat("- GSE256490 (scRNA reference) is large (606,380 cells, 117 samples) - make\n")
cat("  sure you have disk space and a stable connection before 01_download_geo.R.\n")

cat("\nPython packages needed separately (pip install --user):\n")
cat("  requests pandas shap scanpy\n")
