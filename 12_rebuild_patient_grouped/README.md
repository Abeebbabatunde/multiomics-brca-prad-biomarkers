# 12_rebuild_patient_grouped: rebuilding the panels with patients kept together

These scripts are for the follow-up study. They reuse the thesis methods (paired moderated-t filter, LASSO, SVM-RFE and Boruta stability selection, consensus of at least two of three methods) and inputs built by `00_define_cohorts_and_inputs.R`.

**Status: run.** These scripts produced the results reported in thesis Step 12.0 (sections 12.0 to 12.6): cohorts and inputs, stability selection, patient-grouped DIABLO, the leakage demonstrations, the comparison with the initial panels, direction and subtype checks, external validation and the clinical and survival analyses. Run order: `00_define_cohorts_and_inputs.R`, `00_check_patient_groups.R`, `00b_audit_pairing.R`, then `01` (or `01_stability_selection_rebuild_FULL.R`), `02`, `02b`, `02c`, `03`, `04`, `05` and `06`. Check the output of the first three before the long runs; test on one layer (`layers_to_run`) before the full run.

| Script | What it does | Main outputs |
|---|---|---|
| `rebuild_functions.R` | Shared selection functions (same as folder 11) | none |
| `00_define_cohorts_and_inputs.R` | Defines the cohorts (breast stage from the clinical file) and builds harmonised inputs for both cancers: common methylation probes, RNA by vst with blind = TRUE, common genes. Original files untouched | `12_rebuild_inputs/*_rebuild.*` |
| `00_check_patient_groups.R` | Checks that each patient has a tumour and a normal sample in every layer | `12_rebuild_inputs/rebuild_00_input_check.csv` |
| `00b_audit_pairing.R` | Read-only audit of the cohort manifests: codes, patient IDs, vials per patient, same patients in every layer | console output |
| `01_stability_selection_rebuild.R` | Repeats the whole selection on random subsets of patients; selection frequency per feature; panel = frequency >= 0.60; shuffled-label control permuted within patient | `<layer>_selection_frequency_rebuild.csv`, `<layer>_rebuilt_panel.csv`, `*_SHUFFLED.csv` |
| `02_breast_diablo_patient_grouped.R` | Breast DIABLO with patient-grouped folds, keepX grid extended downward, shuffled-label control | `DIABLO_BRCA_patient_grouped_cv.csv`, `DIABLO_BRCA_selected_features.csv` |
| `03_compare_old_and_rebuilt_panels.R` | Compares the locked and rebuilt panels | `rebuild_03_old_vs_rebuilt.csv` |

Put `rebuild_functions.R` in the project folder, or put this whole folder inside it as `12_rebuild_patient_grouped` (script `01` looks in both places). Seed 2026.

Not included: the external re-test of the rebuilt panels. It needs the Step 8.6 external-validation script (GEO downloads and platform mapping) re-run with the rebuilt panel in place of the Step 7 consensus panel, so it should be done once the rebuilt panels exist.


**02b_leakage_demonstration.R** (run after 02): same breast DIABLO model with patient-grouped vs sample-level folds, for keepX 1, 5, 20, 50; writes DIABLO_BRCA_leakage_demonstration.csv.

**02c_selection_leakage_demonstration.R** (run after 02b): same breast DIABLO model, patient-grouped folds throughout, with the label-using feature filter done inside each training fold vs once on all samples; keepX 5 and 20, with shuffled-label controls; writes DIABLO_BRCA_selection_leakage_demonstration.csv.

**04_biological_and_subtype_check.R** (run after 01): paired tumour-minus-normal effect size, direction and consistency of every rebuilt feature, CpG annotation, and PAM50-subtype consistency for breast; writes rebuild_04_effects.csv, rebuild_04_breast_subtype.csv, rebuild_04_cpg_annotation.csv. P-values are descriptive only (features were selected for this contrast).

**05_external_validation_rebuilt_panels.R** (run after 01; needs the GEO series matrices cached by the thesis run): trains a random forest on the rebuild cohort with the rebuilt panels and scores GSE60185 (breast methylation), GSE54002 (breast RNA) and GSE29079 (prostate RNA); writes 12_rebuild_external/*.csv incl. a comparison with the locked panels.

**06_clinical_rebuilt_panels.R** (run after 01): the thesis clinical association (9.1) and survival (9.2) scripts, unchanged except for using the rebuilt panels; writes 12_rebuild_clinical/. **05** now also runs a sensitivity analysis matching CCN4 as WISP1 on the GPL570 array.
