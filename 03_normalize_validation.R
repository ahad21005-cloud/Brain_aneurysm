# =============================================================================
# 03_normalize_validation.R - every validation / blood cohort is normalised on
# its OWN platform and NEVER merged with another cohort or the discovery set.
# Steps: log2 if the 99th percentile > 100 -> quantile normalisation (limma)
# -> probe -> gene symbol (highest-mean probe) -> save.
# GSE75436 additionally stores a verified pair id.
# =============================================================================
suppressPackageStartupMessages({ library(Biobase); library(limma) })
source("scripts/utils/label_utils.R"); source("scripts/utils/gene_mapping.R")
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

normalize_and_map <- function(id, mapping_fn, label_fn, agilent = FALSE, pair_fn = NULL) {
  eset <- readRDS(sprintf("data/raw/%s_eset.rds", id)); expr <- exprs(eset)
  if (quantile(expr, 0.99, na.rm = TRUE) > 100) {
    cat(sprintf("[%s] raw-scale data detected -> log2\n", id)); expr[expr <= 0] <- NA; expr <- log2(expr)
  }
  expr <- normalizeBetweenArrays(expr, method = "quantile")
  map  <- if (agilent) map_gpl4133_agilent(eset) else mapping_fn(rownames(expr))
  expr_gene <- collapse_to_gene_symbol(expr, map)
  status <- label_fn(pData(eset)); check_labels(status, id)
  out <- list(expr = expr_gene, status = status, sample_id = colnames(expr_gene))
  if (!is.null(pair_fn)) out$pair_id <- pair_fn(pData(eset), status)
  saveRDS(out, sprintf("data/processed/%s_normalized.rds", id))
  cat(sprintf("[%s] saved %d genes x %d samples\n\n", id, nrow(expr_gene), ncol(expr_gene)))
}
normalize_and_map("GSE75436", map_gpl570,   label_gse75436, pair_fn = get_pair_id_gse75436)
normalize_and_map("GSE13353", map_gpl570,   label_gse13353)
normalize_and_map("GSE15629", map_gpl6244,  label_gse15629)
normalize_and_map("GSE54083", NULL,         label_gse54083, agilent = TRUE)
normalize_and_map("GSE26969", map_gpl570,   label_gse26969)
normalize_and_map("GSE36791", map_gpl10558, label_gse36791)
cat("Done. Compare label counts with Table 1 of the paper.\n")
