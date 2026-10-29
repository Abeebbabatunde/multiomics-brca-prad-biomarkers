# PRAD proteomics (E-PROT-21), Steps 2.7 and 3.7

| File | Status |
|---|---|
| `02_7_03_7_PRAD_proteomics_AS_RUN.R` | Matches the files on disk: QC, PCA, .rds, .tsv and DEP tables in `03_Results/Proteomics/`, plus the two `write.csv` lines run afterwards. |
| `02_7_03_7_PRAD_proteomics_consolidated_paths_NOT_CONFIRMED_RUN.R` | Same analysis code (only the output paths and the two `write.csv` lines differ). `02_ProcessedData/Prostate/` holds only the .csv, so this version did not produce the other files. |
| `archive/02_7_03_7_template_v1_never_ran.R` | Placeholder paths (`data/proteomics/prostate`); failed with "cannot open file". |

- Input: Expression Atlas export of E-PROT-21 (header "Selected columns: 36 (all)", timestamp 15 September 2026, 10:42:17).
- Raw matrix read with `skip = 4`: 8,483 rows by 38 columns (Gene ID, Gene Name and 36 samples). After filtering: 5,412 proteins (the two intermediate counts were not recorded).
- Design: seed 2026; 28 neoplasm versus 8 adjacent normal samples; limma `~0+group` with `eBayes(trend = TRUE)`; |log2FC| >= 1 and FDR < 0.05 give 4 proteins.
- Individuals 1 to 8 each have a neoplasm and an adjacent-normal sample (from the column headers), but the model does not include the pairing.
- An earlier dead end: `01_RawData/Prostate/prad_prot_meta` (110 TCGA-PRAD RPPA samples, one group level, "EarlyProstateCancer") has no normal samples.
- Files in the raw folder not used by this script: `Human__TCGA_PRAD__MDA__RPPA__..._Firehose_RPPA.cct` and `prad_protein_early.csv`.
