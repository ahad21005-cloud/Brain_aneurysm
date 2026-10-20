# =============================================================================
# 11_figures.R - basic publication figures from the saved result tables.
# (Journal-styled versions in the manuscript were redrawn from the same numbers.)
# Outputs -> results/figures/
# =============================================================================
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
eff <- read.csv("results/tables/signature_score_cohort_effects.csv")
pl  <- read.csv("results/tables/signature_score_pooled.csv")

forest <- function(d, pool, file, title) {
  n <- nrow(d) + 1; png(file, width = 1800, height = 300 + 220 * n, res = 300)
  par(mar = c(4, 9, 3, 1))
  xl <- range(c(d$lo, d$hi, pool$lo, pool$hi)); plot(NA, xlim = xl, ylim = c(0.5, n + 0.5), yaxt = "n", ylab = "",
       xlab = "Hedges' g (lower score  <-  0  ->  higher score)", main = title)
  abline(v = 0, lty = 2, col = "grey40")
  y <- rev(seq_len(nrow(d))) + 1
  segments(d$lo, y, d$hi, y); points(d$g, y, pch = 15)
  segments(pool$lo, 1, pool$hi, 1, lwd = 3, col = "firebrick"); points(pool$g, 1, pch = 18, cex = 2, col = "firebrick")
  axis(2, at = c(y, 1), labels = c(d$cohort, "Pooled"), las = 1); dev.off()
}
p <- eff[eff$analysis == "primary", ]; forest(p, pl[1, ], "results/figures/forest_primary_rupture.png", "ME2 score: ruptured vs unruptured")
a <- eff[eff$analysis == "aneurysm_vs_control", ]; forest(a, pl[2, ], "results/figures/forest_aneurysm_vs_control.png", "ME2 score: aneurysm vs control")

# classifier null (written by 05b)
f <- "results/tables/permutation_null_corrected_v2.csv"; s <- "results/tables/permutation_summary_corrected_v2.csv"
if (file.exists(f) && file.exists(s)) {
  nul <- read.csv(f)$permuted_AUC; sm <- read.csv(s)
  png("results/figures/permutation_null.png", width = 1800, height = 1350, res = 300)
  hist(nul, breaks = 30, col = "grey85", border = "grey40", xlab = "Mean fold AUC", main = "Label-permutation null",
       xlim = range(c(nul, sm$observed_AUC))); abline(v = sm$observed_AUC, col = "red", lwd = 2)
  legend("topleft", bty = "n", col = "red", lwd = 2, legend = sprintf("Observed %.3f (p = %.3f)", sm$observed_AUC, sm$p_value)); dev.off()
}
cat("Figures written to results/figures/\n")
