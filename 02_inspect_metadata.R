# =============================================================================
# 02_inspect_metadata.R
# THE MOST IMPORTANT MANUAL STEP IN THE PIPELINE - do not skip or rush this.
#
# GSE122897 (discovery) matters most here: it's the ONLY cohort carrying all
# three classes at once, so a labeling mistake here corrupts everything
# downstream. Get a second team member to independently check your reading
# of this cohort's labels before moving to Stage 3.
# =============================================================================

library(Biobase)

inspect <- function(rds_path) {
  eset <- readRDS(rds_path)
  cat("\n=====================================================\n")
  cat(basename(rds_path), "\n")
  cat("=====================================================\n")
  pd <- pData(eset)
  cat("Columns available:\n"); print(colnames(pd))
  cat("\n-- title --\n"); print(pd$title)
  char_cols <- grep("characteristics_ch1", colnames(pd), value = TRUE)
  for (cc in char_cols) {
    cat("\n--", cc, "--\n")
    print(table(pd[[cc]], useNA = "always"))
  }
  invisible(pd)
}

cat("\n############# DISCOVERY (highest priority - check twice) #############\n")
pd_122897 <- inspect("data/raw/GSE122897_eset.rds")

cat("\n############# VALIDATION COHORTS #############\n")
pd_75436 <- inspect("data/raw/GSE75436_eset.rds")
pd_13353 <- inspect("data/raw/GSE13353_eset.rds")
pd_15629 <- inspect("data/raw/GSE15629_eset.rds")
pd_54083 <- inspect("data/raw/GSE54083_eset.rds")
pd_26969 <- inspect("data/raw/GSE26969_eset.rds")
pd_6551  <- inspect("data/raw/GSE6551_eset.rds")

cat("\n############# BLOOD #############\n")
pd_36791 <- inspect("data/raw/GSE36791_eset.rds")

cat("\n=====================================================\n")
cat("NEXT STEP: fill in scripts/utils/label_utils.R for every cohort above\n")
cat("based on what was just printed. GSE122897 gets a second-person check.\n")
cat("=====================================================\n")

saveRDS(list(GSE122897 = pd_122897, GSE75436 = pd_75436, GSE13353 = pd_13353,
             GSE15629 = pd_15629, GSE54083 = pd_54083, GSE26969 = pd_26969,
             GSE6551 = pd_6551, GSE36791 = pd_36791),
        "data/metadata/raw_pdata_snapshot.rds")
