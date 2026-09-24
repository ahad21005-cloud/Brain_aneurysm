writeLines('# =============================================================================
# 03_normalize_validation.R  (Role 2 owns this)
# QC + normalize every validation and blood cohort SEPARATELY, each on its
# own platform. None of these are ever merged with each other or with the
# discovery cohort.
#
# GSE6551 intentionally excluded - see label_utils.R comment.
# =============================================================================

library(Biobase)
library(limma)
source("scripts/utils/label_utils.R")
source("scripts/utils/gene_mapping.R")

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

normalize_and_map <- function(eset, mapping_fn, label_fn, cohort_name, is_agilent = FALSE) {
  expr <- exprs(eset)

  qx <- quantile(expr, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE)
  if (qx[6] > 100) {
    cat(sprintf("[%s] Non-log data detected (max=%.1f) - log2 transforming.\\n", cohort_name, qx[6]))
    expr[expr <= 0] <- NA
    expr <- log2(expr)
  }

  expr <- normalizeBetweenArrays(expr, method = "quantile")

  if (is_agilent) {
    probe_to_symbol <- map_gpl4133_agilent(eset)
  } else {
    probe_to_symbol <- mapping_fn(rownames(expr))
  }
  expr_gene <- collapse_to_gene_symbol(expr, probe_to_symbol)

  status <- label_fn(pData(eset))
  check_labels(status, cohort_name)

  out <- list(expr = expr_gene, status = status, sample_id = colnames(expr_gene))
  saveRDS(out, sprintf("data/processed/%s_normalized.rds", cohort_name))
  cat(sprintf("[%s] Saved: %d genes x %d samples\\n\\n", cohort_name, nrow(expr_gene), ncol(expr_gene)))
  out
}

cat(">> GSE75436 (GPL570)\\n")
normalize_and_map(readRDS("data/raw/GSE75436_eset.rds"), map_gpl570, label_gse75436, "GSE75436")

cat(">> GSE13353 (GPL570)\\n")
normalize_and_map(readRDS("data/raw/GSE13353_eset.rds"), map_gpl570, label_gse13353, "GSE13353")

cat(">> GSE15629 (GPL6244)\\n")
normalize_and_map(readRDS("data/raw/GSE15629_eset.rds"), map_gpl6244, label_gse15629, "GSE15629")

cat(">> GSE54083 (GPL4133 Agilent)\\n")
normalize_and_map(readRDS("data/raw/GSE54083_eset.rds"), NULL, label_gse54083, "GSE54083", is_agilent = TRUE)

cat(">> GSE26969 (GPL570)\\n")
normalize_and_map(readRDS("data/raw/GSE26969_eset.rds"), map_gpl570, label_gse26969, "GSE26969")

cat(">> GSE36791 (GPL10558, blood)\\n")
normalize_and_map(readRDS("data/raw/GSE36791_eset.rds"), map_gpl10558, label_gse36791, "GSE36791")

cat("=== Stage 3 (validation) complete. 6 cohorts normalized independently. ===\\n")
', con = "scripts/03_normalize_validation.R")