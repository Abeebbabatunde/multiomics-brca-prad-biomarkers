# Step 1.4: LinkedOmics/Broad iTRAQ proteome (abandoned protein track)

This track is not used in any reported result. BRCA protein in the thesis is RPPA (see `../01_4_BRCA_RPPA_download_and_long_format.R`).

File: `Human__TCGA_BRCA__BI__Proteome__QExact__01_28_2016__BI__Gene__CDAP_iTRAQ_UnsharedLogRatio_r2.cct` (tab-separated; `attrib_name` plus TCGA.xx.xxxx columns; log-ratios; NA allowed). Its header has 105 TCGA sample columns.

| Script | Content |
|---|---|
| `01_4_BRCA_iTRAQ_filter_cct_exact_filename.R` | Loads the .cct file, replaces `.` with `-`, matches the first 12 characters to `clinical_brca.csv$submitter_id`, and writes `01_RawData/Breast/brca_protein_early.csv` and `brca_prot_groups.rds`. It uses the exact file name and is probably the last version; no console output was recorded for it. |
| `01_4_BRCA_iTRAQ_validation_checks.R` | Validation checks (needs `brca_protein_early` in memory). |
| `archive/01_4_BRCA_iTRAQ_filter_all_files_index11.R` | Same logic, selecting the file as `list.files()[11]` (earlier workaround). It also deletes `brca_protein_raw.csv`, `brca_protein_early.csv` and `brca_meth_matrix.csv` before saving. |
| `archive/01_4_early_stage_filter_approaches_1_to_4.R` | Earliest attempt (approach 4 would have kept all samples). |
| `archive/console_*.txt` | Consoles; one shows a wrong file (methylation probes, 461,999 by 0) saved as the protein file, which the later scripts corrected. |

## Consequences
1. The file has no normal samples (all columns are tumour), which is why limma failed with "Coefficients not estimable".
2. `01_RawData/Breast/brca_protein_raw.csv` was deleted in this track. The file now carrying that name is the long-format RPPA table from Step 1.4.
3. The sample count after the clinical match was not displayed; the 76 patients in the manifest are the only count.
