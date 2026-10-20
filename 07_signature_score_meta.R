# =============================================================================
# 07_signature_score_meta.R - PRE-SPECIFIED platform-robust analysis.
#   score = for each cohort INDEPENDENTLY: z-score every ME2 hub gene across that
#   cohort's samples (min 100 hub genes present), then mean z per sample.
#   (Within-cohort scaling removes platform scale differences by construction.)
#   Effect size: Hedges' g (unpaired; paired for GSE75436). Per-cohort Wilcoxon
#   (signed-rank if paired). Pooling: DerSimonian-Laird random effects with
#   t-distribution CI (k-1 df), I^2.
# Analyses
#   PRIMARY   ruptured vs unruptured: GSE13353, GSE15629, GSE54083  (random effects)
#   SECONDARY aneurysm vs control   : GSE15629, GSE54083, GSE75436 (paired), GSE26969
#             - all four, random effects (descriptive)
#             - confirmatory fixed-effect: GSE75436 + GSE26969
#             - sensitivity: drop GSE26969 (random effects, k = 3)
#   EXPLORATORY blood SAH vs control: GSE36791
# Paper: primary pooled g = -0.38 (-4.11 to 3.35), I2 = 0.86; aneurysm-vs-control
# g = -1.28 (-2.06 to -0.50), p = 0.014, I2 = 0; confirmatory -1.21; sensitivity -1.40;
# blood g = -0.62 (-1.18 to -0.06), Wilcoxon p = 0.008.
# =============================================================================
source("scripts/utils/meta_utils.R")
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
hubs <- read.csv("results/tables/hub_genes.csv")$gene
MIN_GENES <- 100

score_cohort <- function(expr) {
  present <- intersect(hubs, rownames(expr))
  if (length(present) < MIN_GENES) stop("only ", length(present), " hub genes present")
  Z <- t(scale(t(expr[present, , drop = FALSE])))      # z across samples, per gene
  Z <- Z[apply(Z, 1, function(x) all(is.finite(x))), , drop = FALSE]
  colMeans(Z)
}
load_scores <- function(id) { d <- readRDS(sprintf("data/processed/%s_normalized.rds", id))
  d$score <- score_cohort(d$expr); d$n_hubs <- length(intersect(hubs, rownames(d$expr))); d }

one_cohort <- function(d, grp1, grp2, label, paired = FALSE) {
  s <- d$score; st <- d$status
  x1 <- s[st %in% grp1]; x2 <- s[st %in% grp2]
  if (paired) {                                        # align by verified pair id
    id1 <- d$pair_id[st %in% grp1]; id2 <- d$pair_id[st %in% grp2]
    stopifnot(setequal(id1, id2), length(id1) == length(unique(id1)))
    x1 <- x1[order(id1)]; x2 <- x2[order(id2)]
    e <- hedges_g_paired(x1, x2); p <- wilcox.test(x1, x2, paired = TRUE, exact = FALSE)$p.value
  } else { e <- hedges_g(x1, x2); p <- wilcox.test(x1, x2, exact = FALSE)$p.value }
  ci <- ci95(e["g"], e["se"])
  data.frame(cohort = label, design = if (paired) "paired" else "unpaired", n1 = e["n1"], n2 = e["n2"],
             hubs_present = d$n_hubs, g = e["g"], se = e["se"], lo = ci[1], hi = ci[2], wilcoxon_p = p,
             row.names = NULL)
}
fmt_pool <- function(p, name) data.frame(pool = name, k = p["k"], g = p["g"], lo = p["lo"], hi = p["hi"],
                                         p = p["p"], I2 = if ("I2" %in% names(p)) p["I2"] else NA, row.names = NULL)

# ---- PRIMARY: ruptured vs unruptured (negative g = lower in rupture) --------
prim <- do.call(rbind, lapply(c("GSE13353", "GSE15629", "GSE54083"), function(id)
  one_cohort(load_scores(id), "ruptured", "unruptured", id)))
prim_pool <- pool_random(prim$g, prim$se)
cat("\n=== PRIMARY: ruptured vs unruptured ===\n"); print(prim, digits = 3); print(round(prim_pool, 3))

# ---- SECONDARY: aneurysm vs control (aneurysm - control; negative = lower in aneurysm)
AN <- c("ruptured", "unruptured", "aneurysm")
sec <- rbind(
  one_cohort(load_scores("GSE15629"), AN, "control", "GSE15629"),
  one_cohort(load_scores("GSE54083"), AN, "control", "GSE54083"),
  one_cohort(load_scores("GSE75436"), AN, "control", "GSE75436", paired = TRUE),
  one_cohort(load_scores("GSE26969"), AN, "control", "GSE26969"))
sec_pool  <- pool_random(sec$g, sec$se)
conf_pool <- pool_fixed(sec$g[sec$cohort %in% c("GSE75436", "GSE26969")], sec$se[sec$cohort %in% c("GSE75436", "GSE26969")])
sens_pool <- pool_random(sec$g[sec$cohort != "GSE26969"], sec$se[sec$cohort != "GSE26969"])
cat("\n=== SECONDARY: aneurysm vs control ===\n"); print(sec, digits = 3)
cat("All four (random):\n"); print(round(sec_pool, 3))
cat("Confirmatory fixed (GSE75436+GSE26969):\n"); print(round(conf_pool, 4))
cat("Sensitivity (drop GSE26969, random):\n"); print(round(sens_pool, 3))

# ---- EXPLORATORY: blood ------------------------------------------------------
bl <- one_cohort(load_scores("GSE36791"), "SAH", "control", "GSE36791")
cat("\n=== EXPLORATORY blood: SAH vs control ===\n"); print(bl, digits = 3)

# ---- save --------------------------------------------------------------------
write.csv(rbind(cbind(analysis = "primary", prim), cbind(analysis = "aneurysm_vs_control", sec),
                cbind(analysis = "blood", bl)), "results/tables/signature_score_cohort_effects.csv", row.names = FALSE)
write.csv(rbind(fmt_pool(prim_pool, "primary RE (k=3)"), fmt_pool(sec_pool, "aneurysm-vs-control RE (k=4)"),
                fmt_pool(conf_pool, "confirmatory FE (GSE75436+GSE26969)"),
                fmt_pool(sens_pool, "sensitivity RE (drop GSE26969)")),
          "results/tables/signature_score_pooled.csv", row.names = FALSE)
cat("\nSaved results/tables/signature_score_cohort_effects.csv and signature_score_pooled.csv\n")
