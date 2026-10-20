# Intracranial aneurysm (IA) transcriptomics: ME2 module, classifier and cross-cohort replication

Code for the MRT Batch A-1 final-year project (Rehman Medical College / Rehman Medical Institute, Peshawar).
All data are public GEO datasets; no patient-identifiable data are used.

## What the analysis shows (final manuscript version)
| Question | Result |
|---|---|
| WGCNA on GSE122897 (58 samples) | 13 modules; ME2 (1,033 genes) status-associated, r = −0.492, p = 8.74e-5; 372 hub genes, down with severity |
| Ruptured-vs-unruptured classifier (LASSO→RF, nested CV, label-blind 5,000-gene pool) | mean fold AUC 0.728 (SD 0.158); permutation null 0.498, p = 0.020 (9/500) |
| Frozen 4-gene panel | LIMCH1, ANKRD29, PLCL1, SHF; best stability 38% (ANKRD29), none ≥ 60% |
| Transfer of frozen model | **Does not transfer** across platforms (non-overlapping raw ranges) |
| ME2 signature score, ruptured vs unruptured (k=3) | pooled g −0.38 (−4.11 to 3.35), I² 0.86 — **no reliable effect** |
| ME2 signature score, aneurysm vs control (k=4) | pooled g −1.28 (−2.06 to −0.50), I² 0 — consistent |
| Blood (GSE36791), exploratory | g −0.62 (−1.18 to −0.06) |
| Gene-level replication | 46/316 genes (14.6%) concordant in all 3 cohorts (chance 12.5%) |
| DGIdb | 94/372 hubs have interactions; none for the 4 panel genes |

Earlier versions of this repo reported a 7-gene panel and AUC 0.708 — that analysis contained a feature-selection leak and an AUC-direction error and is **superseded**.

## Run order (from the project root)
| Script | Purpose |
|---|---|
| `scripts/00_setup.R` | install packages |
| `01_download_geo.R` | download GEO series + supplementary files |
| `02_inspect_metadata.R` | dump sample metadata to `data/metadata/` for label checking |
| `utils/label_utils.R`, `gene_mapping.R`, `classifier_utils.R`, `meta_utils.R` | helpers |
| `03_normalize_discovery.R` | GSE122897: DESeq2 VST; writes sample audit table |
| `03_normalize_validation.R` | each array cohort normalised separately (never merged) |
| `04_wgcna_modules.R` | WGCNA (power 10, signed, minModuleSize 30, merge 0.25) → hub genes |
| `05_ml_classifier.R` | corrected nested CV + frozen panel/model |
| `05b_permutation_test.R` | 500-permutation null, p = (k+1)/(n+1) |
| `05c_stability_selection.R` | 100-bootstrap selection frequencies |
| `06_cross_platform_transfer.R` | frozen model on validation cohorts + range diagnostic |
| `07_signature_score_meta.R` | within-cohort ME2 score, Hedges' g, random/fixed-effect pooling |
| `07b_cohort_independence.R` | shared-sample check |
| `08_gene_level_replication.R` | t-statistic replication + direction concordance |
| `09_scRNA_deconvolution.R` | Seurat reference (GSE256490) + NNLS deconvolution + panel-gene correlations |
| `10_drug_screen_dgidb.py` | DGIdb v5 GraphQL lookup |
| `11_figures.R` | forest plots / null histogram |

`run_all.sh` chains them. Set seeds are fixed (42); record `sessionInfo()` (written by 05/05b).

## Before you trust a run — things to verify
1. **Sample labels** (`utils/label_utils.R`) are keyword rules on GEO metadata. Compare printed label counts with Table 1 of the manuscript and edit the rules if they differ.
2. **GSE122897: 60 GEO samples vs 58 analysed.** `results/tables/GSE122897_sample_audit.csv` lists every sample with its label and reason; use it to state the exact exclusion reason in the manuscript.
3. **GSE75436 pairing**: `get_pair_id_gse75436()` stops unless every pair has exactly one aneurysm and one control sample — adjust its regex after viewing `pdata_GSE75436.csv`.
4. Probe→gene collapsing keeps the highest-mean probe per symbol; state this in Methods.
5. Exact AUC/p/panel can move slightly with R and package versions; the manuscript numbers come from the runs documented in `results/tables/`.
6. `10_drug_screen_dgidb.py` was written against the public DGIdb GraphQL API but not executed by the author of the script's template; check it runs.
