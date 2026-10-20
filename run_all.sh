#!/usr/bin/env bash
# Run the whole pipeline from the project root. Stops at the first error.
set -e
Rscript scripts/00_setup.R
Rscript scripts/01_download_geo.R
Rscript scripts/02_inspect_metadata.R     # then CHECK data/metadata/pdata_*.csv vs label_utils.R
Rscript scripts/03_normalize_discovery.R
Rscript scripts/03_normalize_validation.R
Rscript scripts/04_wgcna_modules.R
Rscript scripts/05_ml_classifier.R
Rscript scripts/05b_permutation_test.R    # slow (500 permutations)
Rscript scripts/05c_stability_selection.R
Rscript scripts/06_cross_platform_transfer.R
Rscript scripts/07_signature_score_meta.R
Rscript scripts/07b_cohort_independence.R
Rscript scripts/08_gene_level_replication.R
Rscript scripts/09_scRNA_deconvolution.R
python scripts/10_drug_screen_dgidb.py
Rscript scripts/11_figures.R
