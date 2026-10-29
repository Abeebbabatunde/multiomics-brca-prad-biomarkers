# Step 1: data downloads

| Script | Content |
|---|---|
| `01_1_BRCA_clinical_AS_RUN.R`, `01_2_BRCA_RNA_AS_RUN.R`, `01_3_BRCA_methylation_AS_RUN.R` | BRCA clinical, RNA and methylation downloads, reconstructed from the code fragments of the session. The run order and the `saveRDS` line for `brca_meth.rds` were not recorded. |
| `01_4_BRCA_RPPA_download_and_long_format.R` | BRCA RPPA: selection of the 30 early-stage pairs, download of the 60 files (GDCdata/TCGA-BRCA/Proteome_Profiling, timestamped 19 June 2026, 12:31) and the long-format table. |
| `01_5_BRCA_metabolomics_exclusion_note.R` | Records the decision to exclude metabolomics. |
| `01_6_PRAD_clinical_AS_RUN.R`, `01_7_PRAD_RNA_AS_RUN.R`, `01_8_PRAD_methylation_AS_RUN.R` | PRAD clinical, RNA and methylation downloads. |

Not included: `00_STARTUP.R` (its writer is in `../00_setup/`) and a PRAD metabolomics download (PRAD metabolomics was excluded).

## Notes
- The RNA and methylation queries use patient barcodes, so tumour, normal and any other sample types of the 804 early-stage patients are all downloaded. The group label in Step 2 treats every sample that is not "Solid Tissue Normal" as tumour.
- The pair selection in `01_4` defines early as Stage I, IA, IB, II, IIA or IIB, keeps all Solid Tissue Normal and Primary Tumor files of early-stage patients, and retains patients with both: 30 patients (33 normal files; 3 belonged to patients without an early-stage primary tumour). One metastatic file (61 rows) was dropped before download. The code contains no filter on vital status or survival.
- `../00_setup/` holds the `00_STARTUP.R` writer and the package and folder checks. The startup file uses outdated step labels (for example "1.4 Prostate clinical data") and `setwd("MyThesis_MultiOmics")`, which fails if the session is already in that folder.
- `archive/` holds failed or superseded attempts (JCI metabolomics, TCGA-RPPA firehose for PRAD without normals, CPTAC quality control). They are not part of the result chain.
- `iTRAQ_track/` holds the abandoned LinkedOmics/Broad iTRAQ protein track (see its README).
