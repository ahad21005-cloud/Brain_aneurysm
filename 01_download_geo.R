# =============================================================================
# 01_download_geo.R
# Downloads every cohort in the FINAL design and caches raw objects to
# data/raw/. Run once per cohort; re-run only if you suspect GEO updated
# a record.
#
# ROLE SPLIT: Role 1 (Discovery Data Lead) only needs the GSE122897 block.
# Role 2 (Validation Data Lead) needs every other block except GSE256490.
# Role 5 (Mechanism Lead) needs the GSE256490 block. Everyone can run this
# script and it will just re-download what's missing - no need to hand-edit
# it per person unless you want to skip sections you don't own.
# =============================================================================

library(GEOquery)
options(timeout = 1200)  # GSE256490 in particular can be slow/large

raw_dir <- "data/raw"
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)

get_and_save <- function(accession, expected_n = NULL, expected_platform = NULL) {
  cat(sprintf(">> Downloading %s...\n", accession))
  gse <- getGEO(accession, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = raw_dir)
  eset <- gse[[1]]
  saveRDS(eset, file.path(raw_dir, paste0(accession, "_eset.rds")))
  n <- ncol(Biobase::exprs(eset))
  plat <- Biobase::annotation(eset)
  ok_n <- is.null(expected_n) || n == expected_n
  ok_p <- is.null(expected_platform) || plat == expected_platform
  flag <- if (ok_n && ok_p) "OK" else "CHECK THIS - mismatch vs expected"
  cat(sprintf("   samples=%-4d platform=%-10s [%s]\n", n, plat, flag))
  invisible(eset)
}

# ---- DISCOVERY (Role 1) -----------------------------------------------------
cat("\n===== DISCOVERY COHORT =====\n")
get_and_save("GSE122897", expected_n = 59, expected_platform = "GPL16791")

# ---- VALIDATION - tissue (Role 2) -------------------------------------------
cat("\n===== VALIDATION COHORTS (tissue) =====\n")
get_and_save("GSE75436",  expected_n = 30, expected_platform = "GPL570")
get_and_save("GSE13353",  expected_n = 19, expected_platform = "GPL570")
get_and_save("GSE15629",  expected_n = 19, expected_platform = "GPL6244")
get_and_save("GSE54083",  expected_n = 23, expected_platform = "GPL4133")
get_and_save("GSE26969",  expected_n = 6,  expected_platform = "GPL570")
# GSE6551 - confirm composition (literature reports 3 UIA + 2 RIA + 5 STA = 10)
# against the actual metadata before trusting it in 02_inspect_metadata.R
get_and_save("GSE6551",   expected_n = 10, expected_platform = "GPL570")

# ---- VALIDATION - blood (Role 2) --------------------------------------------
cat("\n===== BLOOD VALIDATION =====\n")
get_and_save("GSE36791",  expected_n = 61, expected_platform = "GPL10558")

# ---- scRNA REFERENCE (Role 5) -----------------------------------------------
cat("\n===== scRNA DECONVOLUTION REFERENCE =====\n")
cat("GSE256490 (Walchli et al. 2024, Nature) - 606,380 cells, 117 samples.\n")
cat("This is large. The series matrix alone won't have the count data - use\n")
cat("supplementary files.\n")
tryCatch({
  getGEOSuppFiles("GSE256490", baseDir = raw_dir, makeDirectory = TRUE)
}, error = function(e) {
  cat("Supplementary download failed or timed out. Check\n")
  cat("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE256490 manually,\n")
  cat("or check the paper's Nature supplementary data links as a fallback.\n")
  print(e)
})

cat("\n=== All downloads attempted. Review the [OK] / [CHECK THIS] flags above\n")
cat("before proceeding to 02_inspect_metadata.R. ===\n")
