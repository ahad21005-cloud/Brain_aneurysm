# =============================================================================
# 07b_cohort_independence.R - check that no sample appears in two cohorts that
# are pooled: (1) shared GSM accessions, (2) near-identical expression profiles
# across cohorts (Pearson r > 0.999 on shared genes) which would indicate the
# same specimen re-deposited. Run after 03_*. Paper: cohorts independent.
# =============================================================================
suppressPackageStartupMessages(library(Biobase))
ids <- c("GSE122897", "GSE13353", "GSE15629", "GSE54083", "GSE75436", "GSE26969", "GSE36791")
gsm <- lapply(ids, function(i) rownames(pData(readRDS(sprintf("data/raw/%s_eset.rds", i))))); names(gsm) <- ids
ov <- 0
for (a in 1:(length(ids) - 1)) for (b in (a + 1):length(ids)) {
  s <- intersect(gsm[[a]], gsm[[b]]); if (length(s)) { ov <- ov + length(s); cat("SHARED GSM", ids[a], ids[b], ":", s, "\n") } }
cat("Shared GSM accessions across cohorts:", ov, "\n")

expr <- lapply(setdiff(ids, "GSE36791"), function(i) {
  f <- if (i == "GSE122897") "discovery" else i; readRDS(sprintf("data/processed/%s_normalized.rds", f))$expr })
names(expr) <- setdiff(ids, "GSE36791"); hits <- 0
for (a in 1:(length(expr) - 1)) for (b in (a + 1):length(expr)) {
  g <- intersect(rownames(expr[[a]]), rownames(expr[[b]]))
  r <- cor(expr[[a]][g, ], expr[[b]][g, ], method = "spearman")
  if (any(r > 0.999)) { hits <- hits + sum(r > 0.999); cat("Near-identical profiles:", names(expr)[a], names(expr)[b], "\n") } }
cat("Near-identical cross-cohort sample pairs (rho > 0.999):", hits, "\n")
