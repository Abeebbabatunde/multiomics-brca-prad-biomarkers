# Pre-submission re-analyses (nested resampling, external confidence intervals, prostate RNA scale, probe audit)

These scripts reproduce the four corrections made to the TCGA-BRCA / TCGA-PRAD multi-omics biomarker analysis.
They were assembled from the code blocks that were run interactively in R 4.6.0 on Windows 11. Edit the project folder
path at the top of each script (`C:/Users/USER/Documents/MyThesis_MultiOmics`) before running.

| Script | What it does | Main outputs |
| --- | --- | --- |
| `00_build_inputs.R` | Builds the four input pairs (matrix + sample sheet) from the processed data. Methylation: beta values clipped to 0.001-0.999, then M-values; probes with any missing value in the paired samples are dropped. Prostate RNA: log2(normalised count + 1) after keeping genes with count >= 10 in >= 5 samples. | `X_methylation_PRAD.rds`, `X_methylation_BRCA.rds`, `X_rna_BRCA.rds`, `PRAD_rna_log2.rds`, sample sheets |
| `01_nested_patient_cv.R` | Fix 1. Nested, patient-grouped resampling (4 outer folds, 10 repeats) for the four layers. All label-using steps are repeated inside each training fold. Shuffled-label control permutes labels within patients (5 repeats). | `nested_*_summary.csv`, `*_out_of_fold_predictions.csv`, `*_selection_frequency.csv` |
| `02_external_auc_ci.R` | Fix 2. DeLong and stratified-bootstrap 95% CIs for external AUCs; exact intervals for sensitivity and specificity at the locked cut-off. | `external_validation_with_CIs.csv` |
| `03a_prostate_diablo_patient_grouped.R` | Fix 3. DIABLO with patient-grouped folds, log2 RNA, keepX grid 2-50, shuffled-label control. | `DIABLO_PRAD_patient_grouped_cv.csv`, `DIABLO_PRAD_selected_features_annotated.csv` |
| `03b_prostate_mofa2_log2.R` | Fix 3. MOFA2 with the log2 RNA. Run in a fresh R session; needs Python with `mofapy2`. | `MOFA_PRAD_log2.hdf5`, `MOFA_PRAD_factors.rds` |
| `04_prostate_methylation_probe_audit.R` | Fix 4. Counts and annotation of the probes without a P value. | Console output (counts for the Methods) |

Run order: `00`, then `01` to `04`. `01` takes roughly 0.5-2 hours per layer; set `layers_to_run` to run one layer at a time.
Input files from earlier pipeline steps are required (Step 5.0 manifests, `02_ProcessedData`, `03_Results`, the Step 8 external
prediction files with columns `true`, `prob_early`, `pred`, and the Step 10 panel files).

## Reference results (seed 2026)

* Nested median AUC (range): breast methylation 0.999 (0.986-1.000), breast RNA 1.000 (0.996-1.000), prostate methylation
  0.997 (0.983-1.000), prostate RNA 0.986 (0.934-1.000). Shuffled labels: 0.478, 0.478, 0.498, 0.498.
* External AUC (95% CI, DeLong): 0.950 (0.915-0.985), 0.932 (0.915-0.948), 0.947 (0.901-0.992), 0.971 (0.943-0.999).
* Prostate RNA, log2 scale: MOFA2 RNA view 50.7% of variance (about 0.3% before the correction); Factor 1 paired P = 0.00032.
* DIABLO, patient-grouped, real vs shuffled error: keepX 20: 0.038 vs 0.438; keepX 50: 0.024 vs 0.479.
* Probes: 485,577 on the array; 64,426 without a P value (all absent from the 400,283-probe matrix); 63,993 DMPs of 421,151 tested;
  the 1,248 DMPs absent from the matrix are all on chrX/Y.

## Notes

* `01` calls `stats::predict()` explicitly because mixOmics (script `03a`) replaces `predict()` in the session.
* `03a`: the keepX grid is run in one pass here. In the reported run, the rows for 30 and 50 features were added afterwards, so
  their decimals can differ slightly; the final model uses `final_keepX <- 20` as reported.
* `03b`: automatic Python installation (basilisk) failed on the Windows machine used, so MOFA2 was run with a manually installed
  Python 3.11 and `pip install mofapy2` (`use_basilisk = FALSE`).
* Selection-frequency files for RNA layers contain a `symbol` column (Ensembl IDs mapped with org.Hs.eg.db).
## Software versions

The full `sessionInfo()` is in `sessionInfo_folder11.txt`. It was printed in a new R session after loading the packages used by these scripts (data.table, limma, glmnet, Boruta, e1071, pROC, mixOmics, MOFA2, minfi, org.Hs.eg.db), not in the session in which the scripts were first run. It records the installed versions.

| Component | Version |
|---|---|
| R | 4.6.0 (2026-04-24 ucrt), Windows 11 x64 (build 22000) |
| Bioconductor | 3.23 |
| data.table | 1.18.2.1 |
| limma | 3.68.2 |
| glmnet | 5.0 |
| Boruta | 10.0.0 |
| e1071 | 1.7.17 |
| pROC | 1.19.0.1 |
| mixOmics | 6.36.0 |
| MOFA2 | 1.22.0 |
| minfi | 1.58.0 |
| IlluminaHumanMethylation450kanno.ilmn12.hg19 | 0.6.1 |
| org.Hs.eg.db | 3.23.1 |
| mofapy2 (Python 3.11) | 0.7.5 |

The mofapy2 and Python 3.11 versions come from the notes on the MOFA2 run (above) and were not confirmed by this printout: the check in the new session started a different temporary Python (3.12.14, without mofapy2).
