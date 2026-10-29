# Step 3: differential analysis

| Script | Content |
|---|---|
| `03_1_BRCA_DEGs_DESeq2_AS_RUN.R` | Unshrunken DESeq2, `complete.cases`, adjusted P < 0.05 and |log2FC| > 1: 12,402 DEGs (8,558 up, 3,844 down). DESeq2 outlier replacement refitted 8,566 genes. |
| `03_2_BRCA_methylation_limma_AS_RUN.R` | The header reads "Step 2.2", which is a label only; this is the differential step. MASK_general and NA probe filter (485,577 to 305,553 probes), unpaired limma, global BH: 192,919 probes at adjusted P < 0.05. |
| `03_2b_BRCA_DMP_filter_and_copy.R` | Adjusted P < 0.05 and |delta beta| >= 0.2: 18,565 DMPs. The copy is written to `03_DifferentialAnalysis/Breast/brca_DMPs.csv`; Step 5 reads `03_Results/Epigenomics/brca_DMPs.csv`. The saved `brca_DMPs.csv` is identical to the filtered copy. |
| `03_3_BRCA_DEPs_paired_wilcoxon_AS_RUN.R` | RPPA, 457 antibodies by 30 pairs, BH, |mean difference| > 0.5: 114 significant. Saved to `03_DifferentialExpression/Breast/brca_DEPs.csv`; later steps may read `03_Results/Proteomics/brca_DEPs.csv`. |
| `03_5_PRAD_DEGs_DESeq2_AS_RUN.R` | PRAD DEGs (see below). |
| `03_6_PRAD_DMPs_limma_AS_RUN.R` | PRAD DMPs (see below). |

`archive/` holds the never-run Phase 3 templates and an earlier PRAD DEG draft; `console_logs/` holds the raw consoles.

## 3.5 PRAD DEGs
- `03_5_PRAD_DEGs_DESeq2_AS_RUN.R` uses the full path layout and matches the files on disk, including `04_Figures/prad_DEG_volcano.png`. The volcano section failed once ("could not find function EnhancedVolcano"); EnhancedVolcano was then installed and loaded, and the PNG exists.
- `archive/03_5_PRAD_DEGs_draft_results_prostate_folder.R` is the earlier draft (its outputs are not on disk). A diagnostic script that read `path/to/your_prostate_data.rds` failed on the placeholder path and never ran.
- From the console: `dds_prad` has 219 samples, design `~group` and size factors already present, and `resultsNames()` was empty, so `DESeq()` was run in this step and the fitted object was saved over `02_ProcessedData/Prostate/dds_prad.rds`. Contrast EarlyProstateCancer versus Normal, `alpha = 0.05`, `lfcShrink(type = "ashr")`. DEG = padj < 0.05 and |log2FC| >= 1 on the shrunken values; `prad_DEGs.csv` has 5,449 rows and the full table 43,076.
- The ashr output has no `stat` column, which is why the PRAD tables have none (the BRCA table does).
- Not recorded: the printed counts of up and down DEGs, the group table, and `00_STARTUP.R`.

## 3.6 PRAD DMPs (writer of `prad_DMPs_full_results.csv` and `prad_DMPs.csv`)
- The script sources `00_Scripts/00_STARTUP.R`, which is not included. A console error ("could not find function lmFit") shows that limma was not loaded when that line ran, so limma must be loaded by the startup file or beforehand.
- Method: beta clipped to [1e-6, 1 - 1e-6] and converted to M-values; limma `~group` (Normal versus EarlyProstateCancer, unpaired), default `eBayes`, 50,000-probe chunks; delta_beta = rowMeans(tumour) - rowMeans(normal) computed without `na.rm`; DMP = adjusted P < 0.05 and |delta_beta| >= 0.1.
- Why 63,993 and 79,679 DMPs: delta_beta is NA for any probe with an NA in any sample (126,865 NAs in the full table, from the Step 4 recalculation), and `dplyr::filter` drops NA rows, so those probes could not enter `prad_DMPs.csv`. The recalculation (`04_pathway/04_2`) filled 62,439 of them.
- Not recorded: the console counts of significant, hypermethylated and hypomethylated DMPs printed by this script.
