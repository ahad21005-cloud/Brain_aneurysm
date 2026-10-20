# =============================================================================
# 00_setup.R - install every package used by the pipeline. Run once.
# Tested R version for the paper: record yours with sessionInfo() (see README).
# =============================================================================
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc <- c("GEOquery", "Biobase", "DESeq2", "limma", "WGCNA", "org.Hs.eg.db",
          "hgu133plus2.db", "hugene10sttranscriptcluster.db", "illuminaHumanv4.db",
          "AnnotationDbi")
cran <- c("glmnet", "randomForest", "pROC", "matrixStats", "tidyverse", "readr",
          "Seurat", "Matrix", "nnls", "httr", "jsonlite", "ggplot2")

BiocManager::install(setdiff(bioc, rownames(installed.packages())), update = FALSE, ask = FALSE)
install.packages(setdiff(cran, rownames(installed.packages())))

# Python (drug screen only): pip install requests pandas
cat("Setup finished. Now run scripts/01_download_geo.R\n")
