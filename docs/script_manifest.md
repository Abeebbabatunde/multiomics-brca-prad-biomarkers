# Script manifest

Provenance and status of each script. Labels follow the repository README ("Script labels").

| Step | Script | Status |
|---|---|---|
| 0 | `00_4_write_STARTUP_script`, `00_checks` | Pasted code, filed unchanged |
| 1.1 to 1.3 | `01_1` BRCA clinical, `01_2` BRCA RNA, `01_3` BRCA methylation | Assembled from code fragments; the `saveRDS` line for `brca_meth.rds` was not recorded |
| 1.4 | `01_4` BRCA RPPA download, 30-pair selection and long format | Assembled from the code of the 19 June 2026 session |
| 1.5 | `01_5` metabolomics exclusion note | Pasted code, filed unchanged |
| 1.6 to 1.8 | `01_6` PRAD clinical, `01_7` PRAD RNA, `01_8` PRAD methylation | As run |
| 0.5 | `00_5_manifest_RPPA_EDITED`, `00_5_representativeness_BRCA_proteomics` | Edited after the first run (protein source changed to RPPA) |
| 0.5 (paired) | `00_5_paired_manifest_RECOVERED_FROM_HISTORY` | Recovered from the RStudio history (20 September 2026); run after `05_0a` |
| 2.1 to 2.6 | `02_1`, `02_2`, `02_3`, `02_5`, `02_6`; `02_BRCA_removed_probe_audit` | As run; the probe audit was reconstructed |
| 2.7 / 3.7 | `02_7_03_7` PRAD proteomics | As run (the consolidated-paths version is not confirmed run) |
| 3.1 to 3.6 | `03_1`, `03_2`, `03_2b`, `03_3`, `03_5`, `03_6` | As run |
| 4.1 to 4.5 | `04_1` to `04_5` | Reconciled against the printed results |
| 5.0 to 5.6 | `05_0a`, `05_0b`, `05_1a`, `05_1b`, `05_2`, `05_3`, `05_4a`, `05_4b`, `05_5` (markers, as run), `05_5b` and `05_5c` (whole-symbol markers, recovered from history), `05_6` (Python environment) | Reconciled against the printed results |
| 6.0a | `06_0a_consensus_features.R` | Final |
| 6.0c | `06_0c_confound_check_CORRECTED.R` | Final (the 6.0b logistic version was retired) |
| 6.1 | `06_1_network_pathway.R` | Final (corrected, namespaced) |
| 6.1b | `06_1b_hallmark_nominal.R` | Final |
| 6.2 / 6.2b | `06_2_figures_BRCA_meth_rna.R`, `06_2b_figures_BRCA3omics_PRAD.R` | Final |
| 6.3 | `06_3_integrated_interpretation.R` | Final (no `categorySize`) |
| 6.4 | `06_4_network_statistics.R` | Final |
| 7.0 | `07_0_targeted_feature_selection_ML.R` | Final (patient-level) |
| 8.0 | `08_0_internal_validation.R` | Final |
| 8.6 | `08_6_external_validation.R` | Final (download with retries) |
| 9.0 to 9.2 | `09_0_local_audit.R`, `09_1_clinical_association_RPPA_EDITED.R`, `09_2_survival.R` | Final; 9.1 edited so that VAV1 uses RPPA (n = 30) |
| 10.1 / 10.3b | `10_1_final_panel_tables_EDITED.R`, `10_3b_rna_fc_annotation_pathway_overlap.R` | 10.1 edited and rerun (n = 30); 10.3b is unaffected (it reads Step 7 outputs and excludes VAV1) |
| 11 | `00_build_inputs`, `01_nested_patient_cv`, `02_external_auc_ci`, `03a_prostate_diablo_patient_grouped`, `03b_prostate_mofa2_log2`, `04_prostate_methylation_probe_audit` | Assembled from the code blocks run interactively after Step 10; filed unchanged. In `03a` the keepX rows for 30 and 50 were added after the reported run, so their decimals can differ slightly. `03b` ran with a manually installed Python 3.11 and mofapy2 (`use_basilisk = FALSE`) |

## Not in the result chain
- Superseded and troubleshooting versions (GDC audit, Hedges g version of 9.1, first Step 8 external script, 8.0b inspection, 7.1 to 7.2C, 6.0b, earlier 6.1 to 6.3) are not part of this library or are filed in `archive/` folders.
- `00_setup/archive/00_confirm_outstanding_items.R` is a read-only diagnostic and changes no data files.

## Not recorded
- The `00_Scripts` folder listing and `.Rhistory` were not available, so the project's script folder was not matched against these scripts. Every script comes from project files, the RStudio history or pasted code.
- The console message for the GSE76938 filter (134 of 136 samples retained) was not kept, so the panel values missing for the 2 dropped samples are not identified.

- The `sessionInfo()` of the folder 11 scripts was taken in a new session after the runs (`sessionInfo_folder11.txt`).
- The rule that removed the 64,426 probes without a P value before testing (PRAD Step 1.0) was not recorded.

## Known inconsistencies
- The PRAD RNA scale differs between the original Step 5 (count scale) and the folder 11 re-run (log2).
- The RNA scale differs across steps: Steps 7, 8 and 10 use DESeq2-normalised counts, and Step 9.1 uses variance-stabilised counts for BRCA.
- `Step5_Integration/*_Final.rds` files are logical vectors rather than data objects.
