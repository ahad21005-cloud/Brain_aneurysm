# =============================================================================
# 08_gene_level_replication.R - secondary analysis (a): do the 372 hub genes show
# the same ruptured-vs-unruptured direction in the independent cohorts?
#  * Welch t per gene (ruptured - unruptured) in discovery and each cohort
#  * Spearman rho between the two t-vectors over shared hub genes
#  * one-sided permutation p (1,000 shuffles of the validation t-vector)
#  * % genes lower in rupture; direction concordance with discovery across cohorts
#    (chance for all three = 0.5^3 = 12.5%). Paper: 46/316 = 14.6% concordant in all 3.
# =============================================================================
set.seed(42)
hubs <- read.csv("results/tables/hub_genes.csv")$gene
welch_t <- function(E, st, genes) {
  a <- E[genes, st == "ruptured", drop = FALSE]; b <- E[genes, st == "unruptured", drop = FALSE]
  vapply(seq_along(genes), function(i) tryCatch(unname(t.test(a[i, ], b[i, ])$statistic), error = function(e) NA_real_), numeric(1))
}
disc <- readRDS("data/processed/discovery_normalized.rds")
gd <- intersect(hubs, rownames(disc$expr)); names(gd) <- gd
t_disc <- setNames(welch_t(disc$expr, disc$status, gd), gd)

cohorts <- c("GSE13353", "GSE15629", "GSE54083"); tv <- list(); summ <- list()
for (nm in cohorts) {
  d <- readRDS(sprintf("data/processed/%s_normalized.rds", nm)); keep <- d$status %in% c("ruptured", "unruptured")
  E <- d$expr[, keep]; st <- d$status[keep]; g <- intersect(gd, rownames(E))
  tc <- setNames(welch_t(E, st, g), g); tv[[nm]] <- tc
  ok <- is.finite(tc) & is.finite(t_disc[g]); x <- t_disc[g][ok]; y <- tc[ok]
  rho <- cor(x, y, method = "spearman")
  perm <- replicate(1000, cor(x, sample(y), method = "spearman"))
  summ[[nm]] <- data.frame(cohort = nm, genes_compared = length(x), spearman_rho = rho,
                           perm_p_one_sided = (sum(perm >= rho) + 1) / (1000 + 1),
                           pct_lower_in_rupture = 100 * mean(y < 0))
}
summ <- do.call(rbind, summ); print(summ, digits = 3)
write.csv(summ, "results/tables/gene_level_replication.csv", row.names = FALSE)

# concordance of direction with discovery
all3 <- Reduce(intersect, lapply(tv, names)); n3 <- length(all3)
conc <- sapply(tv, function(tc) sign(tc[all3]) == sign(t_disc[all3]))
nconc <- rowSums(conc)
cat(sprintf("\nGenes tested in all 3 cohorts: %d\n", n3)); print(table(factor(nconc, levels = 0:3)))
cat(sprintf("Concordant in all 3: %d (%.1f%%); chance = 12.5%%\n", sum(nconc == 3), 100 * mean(nconc == 3)))
ct <- data.frame(gene = all3, n_concordant = nconc, t_disc = t_disc[all3], sapply(tv, function(x) x[all3]), row.names = NULL)
ct <- ct[order(-ct$n_concordant), ]
write.csv(ct, "results/tables/gene_level_concordance.csv", row.names = FALSE)
cat("\nTop concordant genes:\n"); print(head(ct, 10), digits = 3)
