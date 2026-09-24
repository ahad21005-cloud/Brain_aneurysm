# =============================================================================
# 07_blood_validation.R  (Role 4 owns this)
# Applies the FROZEN panel + classifier to GSE36791 (peripheral blood,
# SAH vs healthy). Framing reminder: blood was drawn AFTER rupture, so this
# tests post-event detectability, not prediction beforehand. State that
# explicitly in the writeup regardless of the result.
# =============================================================================

library(caret)
library(pROC)

final_model <- readRDS("results/model/final_classifier.rds")
selected_genes <- read.csv("results/tables/final_gene_panel.csv")$gene
blood <- readRDS("data/processed/GSE36791_normalized.rds")

common <- intersect(selected_genes, rownames(blood$expr))
cat(sprintf("Panel genes present in blood cohort: %d / %d\n", length(common), length(selected_genes)))
if (length(common) < length(selected_genes)) {
  cat("Missing from blood platform (Illumina): ", setdiff(selected_genes, common), "\n")
}

X_blood <- t(blood$expr[common, , drop = FALSE])
missing <- setdiff(selected_genes, common)
if (length(missing)) {
  pad <- matrix(0, nrow = nrow(X_blood), ncol = length(missing), dimnames = list(rownames(X_blood), missing))
  X_blood <- cbind(X_blood, pad)[, selected_genes, drop = FALSE]
}

y_blood <- factor(ifelse(blood$status == "SAH", "ruptured", "not_ruptured"), levels = c("not_ruptured", "ruptured"))

pred_prob <- predict(final_model, X_blood, type = "prob")[, "ruptured"]
roc_obj <- roc(y_blood, pred_prob, quiet = TRUE, levels = rev(levels(y_blood)))
auc_val <- as.numeric(auc(roc_obj))
ci_val <- ci.auc(roc_obj)

cat(sprintf("\nBlood validation AUC (SAH vs healthy): %.3f\n", auc_val))
cat(sprintf("95%% CI: %.3f - %.3f\n", ci_val[1], ci_val[3]))

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
png("results/figures/roc_blood_validation.png", width = 500, height = 500)
plot(roc_obj, main = sprintf("Blood validation (GSE36791)\nAUC = %.3f", auc_val))
dev.off()

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(data.frame(cohort = "GSE36791_blood", AUC = auc_val, CI_lower = unname(ci_val[1]),
                     CI_upper = unname(ci_val[3])),
          "results/tables/blood_validation_auc.csv", row.names = FALSE)

cat("\nSaved: results/tables/blood_validation_auc.csv\n")
cat("Saved: results/figures/roc_blood_validation.png\n")
cat("\nIf the CI crosses 0.5, report this as null/inconclusive rather than\n")
cat("overstating a non-significant point estimate.\n")
cat("\n=== Stage 7 complete. ===\n")