# =============================================================================
# 02_inspect_metadata.R - dump the sample table (pData) of every cohort to CSV
# so you can CHECK the group labels by eye before trusting label_utils.R.
# Output: data/metadata/pdata_<GSE>.csv  and  data/metadata/raw_pdata_snapshot.rds
# =============================================================================
suppressPackageStartupMessages(library(Biobase))
dir.create("data/metadata", showWarnings = FALSE, recursive = TRUE)
ids <- c("GSE122897", "GSE13353", "GSE15629", "GSE54083", "GSE75436", "GSE26969", "GSE36791")
snap <- list()
for (id in ids) {
  f <- file.path("data/raw", paste0(id, "_eset.rds"))
  if (!file.exists(f)) { cat("Missing:", f, "\n"); next }
  pd <- pData(readRDS(f)); snap[[id]] <- pd
  write.csv(pd, file.path("data/metadata", paste0("pdata_", id, ".csv")))
  cat(sprintf("\n===== %s: %d samples =====\n", id, nrow(pd)))
  chr <- pd[, grepl("title|characteristics|source_name", colnames(pd)), drop = FALSE]
  print(head(chr, 4))
}
saveRDS(snap, "data/metadata/raw_pdata_snapshot.rds")
cat("\nOpen the pdata_*.csv files and confirm the labels match the GEO pages.\n")
