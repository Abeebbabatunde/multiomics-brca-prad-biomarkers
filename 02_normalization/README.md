# Step 2: normalisation

| Script | Content | Counts |
|---|---|---|
| `02_1_BRCA_RNA_AS_RUN.R` | Genes with reads in at least 10% of samples; DESeq2 size factors; `vst(blind = FALSE)`; PCA | see the console log |
| `02_2_BRCA_methylation_quantile_AS_RUN.R` | Probes with less than 20% missing values; chrX and chrY removed; quantile normalisation | 393,747 probes |
| `02_3_BRCA_RPPA_QC_median_impute_AS_RUN.R` | Targets with more than 20% missing values dropped; per-target median imputation | 457 antibodies, 30 pairs; writes `02_Normalized/Breast/brca_protein_normalized.csv` (long format) |
| `02_5_PRAD_RNA_AS_RUN.R` | As 2.1 without vst; size-factor-normalised counts | not recorded |
| `02_6_PRAD_methylation_quantile_AS_RUN.R` | As 2.2 | 400,283 probes |
| PRAD proteomics (2.7) | See `../02_03_proteomics_PRAD/` | |

Metabolomics steps 2.4 and 2.8 are not included (metabolomics was excluded).

## Notes
- `archive/02_3_BRCA_iTRAQ_normalisation_FIX_not_used_in_Step5.R` normalises the separate LinkedOmics iTRAQ proteome (`02_ProcessedData/Breast/brca_protein_normalized.csv`; 50% missing cut-off; minimum/2 imputation). The BRCA protein layer used in Step 5 is the RPPA file above, not this one.
- RNA and methylation group labels treat every sample that is not "Solid Tissue Normal" as tumour. The sample-type composition of the 905 BRCA RNA columns and 219 PRAD RNA columns was not displayed.
