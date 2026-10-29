# Step 0.5: patient manifest and representativeness

| File | Status |
|---|---|
| `00_5_manifest_RPPA_EDITED.R` | Builds the patient manifest and attrition report. BRCA protein is the RPPA layer (`02_Normalized/Breast/brca_protein_normalized.csv`; patient IDs from `sample_barcode`). PRAD `complete_clinical` uses `ajcc_pathologic_stage`, which is empty for PRAD, so PRAD eligible = 0 in the printed attrition. |
| `00_5_manifest_gleason_fix_OUTPUT_NOT_SEEN.R` | Same script with PRAD `complete_clinical` based on `gleason_score`. This one line is the only difference. Its output was not recorded, and no run of it appears in the RStudio history from 20 September 2026 onward. |
| `00_5_representativeness_BRCA_proteomics.R` | Compares BRCA patients with and without RPPA data (chi-square, Fisher's exact test, tableone). |
| `00_5_paired_manifest_RECOVERED_FROM_HISTORY.R` | Code that wrote `sample_manifest_paired_normal_early.csv` and `attrition_report_paired_design.csv` on 20 September 2026 at 16:32:48 (copied verbatim from the RStudio history). It merges the two BRCA cohorts from `05_0a` (16 and 51 pairs) with stage, vital status and follow-up from the harmonised manifest. Its attrition numbers (80, 64, 28, 16, 51) are typed in the script. Run it after `05_0a`; Step 7 reads its output. |
| `archive/00_5_manifest_FIRST_RUN_iTRAQ_protein.R` | First-run version, which read the leftover iTRAQ protein file (76 patients, complete omics 31). |
| `archive/00_5_template_v1_never_ran.R` | Template with placeholder paths; it failed with "cannot open 01_epigenomics/dmp_results.rds". |

## Printed attrition (first run, iTRAQ protein)
BRCA 803 / epigenomics 570 / transcriptomics 801 / protein 76 / complete omics 31 / complete clinical 803; PRAD 188 / 188 / 187 / none / 187 / 0.
With the RPPA edit the BRCA protein count is 30 and complete omics 25 (reported in the thesis).

## Notes
- The first-run manifest read `02_ProcessedData/Breast/brca_protein_normalized.csv` (9,406 rows by 76 patient columns). Step 5 (`05_0a`, `05_1a`) reads `02_Normalized/Breast/brca_protein_normalized.csv`, a different file with a `sample_barcode` column.
- No supplied script reads the outputs of `00_5_manifest_RPPA_EDITED.R` (`sample_manifest_harmonized.csv`, `attrition_report.csv`). Step 7 reads `sample_manifest_paired_normal_early.csv`, written by the recovered script above.
- The PRAD cohort is not defined from the manifest. `05_0b` takes Early as Gleason score of 7 or below directly from `prad_rna_meta.rds` (identical to the RStudio history of 22 September 2026, 02:17).
- `05_0a` is identical to the cohort block in the RStudio history of 20 September 2026, 16:18, which also wrote `cohort_check.txt`.
- The history blocks at 16:05 to 16:12 and 16:28 on 20 September (`preflight.txt`, `harmonization_check.txt`) are read-only diagnostics and are not part of the result chain.
