# Results

Curated figures used in the thesis. Raw and intermediate outputs are excluded (see `.gitignore`).
Figure numbers refer to the thesis Results chapter. Files marked "raster" are image copies; the original
PNG/PDF files remain in the local project folders (not tracked).

## figures/

| File | Thesis figure | Content |
|---|---|---|
| PCA_Breast_RNA.pdf | 2.1 | PCA of BRCA RNA-seq samples |
| Volcano_Breast_DEGs.pdf | 2.2 | BRCA DEG volcano plot |
| prad_DEG_volcano.jpg | 2.3 | PRAD DEG volcano plot (raster) |
| prad_DEP_volcano_plot.pdf | 3.1 | PRAD DEP volcano plot |
| prad_DEP_heatmap.pdf | 3.2 | PRAD DEP heatmap |
| prad_protein_PCA.pdf | 3.3 | PCA of normalised PRAD proteomics data |
| prad_protein_QC_boxplots.pdf | 3.4 | PRAD protein abundance before and after quantile normalisation |

## figures/step5/ (Step 5.0 integration; BRCA = 51 methylation-RNA pairs, PRAD = 17 pairs)

| File | Thesis figure |
|---|---|
| BRCA_upset.jpg | 5.1 |
| PRAD_upset.jpg | 5.2 |
| BRCA_DIABLO_indiv.jpg | 5.3 |
| BRCA_DIABLO_loadings.jpg | 5.4 |
| BRCA_DIABLO_circos.jpg | 5.5 |
| PRAD_DIABLO_indiv.jpg | 5.6 |
| PRAD_DIABLO_loadings.jpg | 5.7 |
| PRAD_DIABLO_circos.jpg | 5.8 |
| BRCA_MOFA_variance.jpg | 5.9 |
| BRCA_MOFA_factors.jpg | 5.10 |
| PRAD_MOFA_variance.jpg | 5.11 |
| PRAD_MOFA_factors.jpg | 5.12 |
| BRCA_venn.jpg, PRAD_venn.jpg | not in thesis (Venn versions of 5.1 and 5.2; set labels are clipped in these images) |

## figures/step6/ (Step 6.0 network and pathway analysis)

| File | Thesis figure |
|---|---|
| BRCA_meth_rna_FULL_network.jpg | 6.1 |
| BRCA_3omics_FULL_network.jpg | 6.2 |
| PRAD_meth_rna_FULL_network.jpg | 6.3 |
| BRCA_meth_rna_Immune_vs_ECM_dotplot.jpg | 6.4 |
| BRCA_meth_rna_cnetplot_GO.jpg | 6.5 |
| BRCA_3omics_nominal_dotplot.jpg | 6.6 |
| PRAD_meth_rna_nominal_dotplot.jpg | 6.7 |
| BRCA_meth_rna_FULL_GO_BP_dotplot.jpg | not in thesis (single-term plot, covered by 6.5) |

## figures/step9/ (Step 9.0 clinical and prognostic validation, BRCA)

| File | Thesis figure | Content |
|---|---|---|
| Figure9_1_composite.png | 9.1 | Panels A to D taken unaltered from pages 1, 2, 3 and 5 of Step9_1_BRCA_significant.pdf |
| Figure9_2_composite.png | 9.2 | Panels A to D taken unaltered from pages 1, 2, 3 and 11 of Step9_2_KM_median_split.pdf |
| Step9_1_BRCA_significant.pdf | repository only | All 27 associations at FDR < 0.05, one per page |
| Step9_2_KM_median_split.pdf | repository only | Median-split Kaplan-Meier curves for all 12 features, one per page |

## tables/

Small result tables (copies of files in the local project folders; large and intermediate tables are excluded).

| File | Source folder | Content |
|---|---|---|
| Step8_FINAL_external_summary.csv | 08_FinalResults/Step8_External_Validation | External AUC per cohort and layer (GEO series, n, features used) |
| Step10_1_BRCA_final_panel.csv | 10_FinalResults/Step10_FinalPanel | BRCA final panel: direction, external, clinical and prognostic evidence per feature |
| Step10_2_PRAD_final_panel.csv | 10_FinalResults/Step10_FinalPanel | PRAD final panel: direction, external, clinical and prognostic evidence per feature |
| Step9_1_clinical_association_results.csv | 09_FinalResults/Step9_Clinical | All 84 clinical association rows (BRCA and PRAD) |
| Step9_1_feature_availability.csv | 09_FinalResults/Step9_Clinical | Number of patients with a value per feature |
| Step10_pathway_overlap.csv | 10_FinalResults/Step10_FinalPanel | Step 6 enrichment files checked for terms shared between cancers |

## Re-analyses after Step 10 (folder `11_pre_submission_reanalyses/`)
No figures were produced by these scripts. Their reported values are given in the thesis (Results, "Nested patient-grouped resampling", "External validation with confidence intervals", "Prostate RNA scale, MOFA2 and DIABLO", "Probes without a test result") and in the reference results of the folder README. Output files are written to the local project folders and are not tracked.
