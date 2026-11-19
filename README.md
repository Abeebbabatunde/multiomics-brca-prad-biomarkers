# MyThesis_MultiOmics: script library

In silico multi-omics biomarker discovery in early-stage breast (TCGA-BRCA) and prostate (TCGA-PRAD) cancer, written in R. Every script here is the code that produced the reported results, filed with a label that says how it reached this library (see "Script labels"). Nothing was written afresh for the library.

## Citation
If you use these scripts, please cite:

Jamiu AB (2026). Multi-omics biomarker discovery in early-stage breast and prostate cancer: analysis scripts. GitHub. https://github.com/Abeebbabatunde/multiomics-brca-prad-biomarkers

`CITATION.cff` holds the same information in machine-readable form (GitHub shows it as "Cite this repository").

## Setup
- R 4.6.0 (2026-04-24), Bioconductor 3.23, Windows 11 (x64).
- Packages: TCGAbiolinks 2.40.0, SummarizedExperiment 1.42.0, minfi 1.58.0, limma 3.68.2, edgeR 4.10.0, DESeq2 1.52.0, MOFA2 1.22.0, mixOmics 6.36.0, glmnet 5.0, randomForest 4.7.1.2, e1071 1.7.17, survival 3.8.6, sesameData 1.30.0. Step 7 also needs kernlab (SVM-RFE). The pre-submission re-analyses (`11_pre_submission_reanalyses/`) also use Boruta, pROC and mofapy2 0.7.5 (Python 3.11); their versions (Boruta 10.0.0, pROC 1.19.0.1) and the full `sessionInfo()` are in `11_pre_submission_reanalyses/`.
- Seeds: 123 (MOFA2, Step 5), 1 (shuffled-label check, Step 5), 42 (Steps 7, 8 and 9.1, and the random forests of the external validation in folder 12), 2026 (all re-analyses in folders 11 and 12). The Step 4 GSEA scripts set no seed.
- Every script starts with `setwd("C:/Users/USER/Documents/MyThesis_MultiOmics")`; change this to your project folder. `00_setup/` writes the `00_STARTUP.R` helper and checks packages and folders.
- The scripts read and write relative to the project folder: `GDCdata`, `01_RawData`, `02_ProcessedData`, `02_Normalized`, `00.5_clinical_metadata_harmonization`, `03_Results`, `05_FinalResults`, `06_FinalResults`, `07_FinalResults` and later result folders. Raw and processed data are not part of this library.

## Run order
The original run order was not recorded. This order follows the step numbers and the dependencies visible in the scripts.

| Step | Folder | What it does |
|---|---|---|
| 0 | `00_setup/` | Startup script writer, package and folder checks |
| 1 | `01_download/` | Clinical, RNA, methylation and RPPA downloads (BRCA); clinical, RNA, methylation (PRAD) |
| 0.5 | `00_5_manifest/` | Patient eligibility manifest, representativeness of the RPPA patients; run once the normalised matrices from Step 2 exist |
| 2 | `02_normalization/`, `02_methylation_probe_audit/`, `02_03_proteomics_PRAD/` | Per-layer filtering and normalisation; PRAD proteomics normalisation and test |
| 3 | `03_differential/` | DEGs, DMPs and DEPs |
| 4 | `04_pathway/` | Individual-omics pathway analysis |
| 5 | `05_integration/` | Paired Normal-versus-Early cohorts, MOFA2 and DIABLO, overlap plots |
| 5 (markers) | `05_integration/05_5b…` and `05_5c…` | Whole-symbol classical-marker scripts, recovered from RStudio history (`collect_needed_2.R`, `collect_needed_3.R`); they also contain read-only checks |
| 0.5 (paired) | `00_5_manifest/00_5_paired_manifest_RECOVERED_FROM_HISTORY.R` | Needs the Step 5.0 cohort files (`05_0a`), so run it after them; Step 7 reads its output |
| 6 | `06_network_pathway/` | Consensus features, confound check, networks, enrichment, figures |
| 7 | `07_feature_selection/` | Targeted feature selection and machine learning |
| 8 | `08_validation/` | Internal and external validation |
| 9 | `09_clinical/` | Clinical association and survival |
| 10 | `10_final_panel/` | Final panel tables and annotation |
| 11 | `11_pre_submission_reanalyses/` | Re-analyses run after Step 10: nested patient-grouped resampling, external AUC confidence intervals, prostate RNA on the log2 scale (MOFA2 and DIABLO), probe audit. Run `00_build_inputs.R` first, then `01` to `04`; the folder README lists inputs and reference results |
| 12 | `12_rebuild_patient_grouped/` | The rebuild with patients kept together: cohorts and harmonised inputs, stability selection with shuffled-label controls, patient-grouped breast DIABLO, leakage demonstrations, comparison with the initial panels, direction and subtype checks, external validation and clinical analyses (thesis Step 12.0). The folder README gives the run order |
| 13 | `13_manuscript_support/` | Scripts that draw the manuscript figures and run the checks made during manuscript preparation |

Each step folder lists its scripts in numeric order (`05_0a`, `05_0b`, `05_1a` and so on).

## Summary of results

| Folders | Main results |
|---|---|
| `01` to `04` (one layer at a time) | **Breast:** 18,565 of 305,553 methylation sites changed at the threshold used (change in methylation of 0.2 or more; 65,320 at 0.1 or more, the threshold used for prostate); 12,402 of 44,607 genes changed in RNA; 114 of 457 proteins changed. **Prostate:** 79,679 of 421,151 methylation sites changed at 0.1 or more (after filling 62,439 missing values of the change in methylation in `04_2`; the original column gives 63,993); 5,449 of 43,076 genes changed in RNA; 4 of 5,412 proteins changed (AMACR higher; MELTF, MYLK and CAPG lower) |
| `05_integration` | MOFA2 and DIABLO in patient pairs: 16 breast patients with three layers, 51 breast and 17 prostate with two layers. In breast, the first factor separated tumour from normal (P = 3.05e-5 with three layers), and methylation and RNA each explained a large share of the variation (62.8% and 58.9% in the 51 pairs) |
| `06_network_pathway` | Only the full breast methylation-and-RNA gene set gave enriched pathways (immune-related: allograft rejection, immunoglobulin-superfamily cell-adhesion signalling and T-cell activation). The breast three-layer set and the prostate set gave none |
| `07` to `10` (panels) | Initial panels of 13 breast and 8 prostate features, with no feature shared between the cancers. The AUC of 1.00 inside TCGA is not a performance estimate. In independent datasets the AUC was 0.932 to 0.971. No breast feature predicted survival; prostate survival could not be assessed (2 deaths) |
| `11_pre_submission_reanalyses` | With every selection step repeated inside each training fold, the held-out AUC stayed high (0.986 to 1.000) but the selected features were unstable: none of the nine breast CpGs was chosen again. External AUCs with 95% intervals: breast methylation 0.950, breast RNA 0.932, prostate methylation 0.947, prostate RNA 0.971. After correcting the RNA scale, RNA explained 50.7% of prostate variation (about 0.3% before) |
| `12_rebuild_patient_grouped` | Rebuilt panels (56 breast and 17 prostate pairs): breast 2 CpGs and 10 genes; prostate 3 genes and no methylation marker. With shuffled labels no feature passed the threshold. Only CAPN11 was in both the original and the rebuilt panels. External AUC 0.949 to 1.000. Breast DIABLO error was 0.006 to 0.035 with patients kept together (0.47 to 0.52 with shuffled labels). No breast feature was prognostic |
| `13_manuscript_support` | Figure scripts, the breast methylation threshold count, and the sex and gender check |

**What the results show.** Tumour and normal tissue differ widely in methylation and RNA in both cancers and in breast protein; in prostate protein only four proteins changed. A very high AUC does not show that a marker panel is reliable: different marker sets gave similar external AUCs, and the markers selected changed when the patients changed. The breast and prostate panels share no marker. The panels discriminate tumour from normal, and no breast feature was prognostic. The cohorts are small (16 to 56 breast pairs and 17 prostate pairs), the panels were derived in TCGA cohorts only, and the rebuilt panels were derived and described in the same patients, so they need testing in independent patients.

## Script labels (end of file name)
The label at the end of a file name records how the script reached this library. It is a record of provenance, not a grade of quality.

| Label | Meaning |
|---|---|
| `AS_RUN` | The script as it was run, matched to the printed results. |
| `AS_PASTED` | Pasted from the session and filed unchanged. |
| `EDITED` | Changed after a first run; only the protein source lines differ (BRCA protein is RPPA throughout). The first-run version is in that folder's `archive/`. |
| `RECOVERED_FROM_HISTORY` | Copied from the RStudio history (20 Sep 2026). |
| `NOT_CONFIRMED_RUN`, `OUTPUT_NOT_SEEN` | Present for completeness; their output was not confirmed. |
| No label | Reconciled against the printed results and the document. |

- `archive/` folders hold failed or superseded attempts and are not part of the result chain. `01_download/iTRAQ_track/` is the abandoned LinkedOmics/Broad iTRAQ protein track; no reported result uses it.
- Folder 11 scripts are filed as supplied and carry no label suffix; they were run interactively in R 4.6.0 (see the folder README for the notes on `03a` and `03b`).

## Data sources
Data are not stored in this repository.
- TCGA-BRCA and TCGA-PRAD methylation (Illumina 450K), expression (STAR counts) and clinical data: NCI Genomic Data Commons, through TCGAbiolinks (BRCA retrieved 9 May 2026, PRAD 18 to 19 May 2026).
- BRCA protein: GDC reverse-phase protein array files (Proteome Profiling, Protein Expression Quantification), downloaded 19 June 2026.
- PRAD protein: E-PROT-21 (Expression Atlas, exported 15 September 2026).
- External validation (GEO): GSE60185 (BRCA methylation), GSE54002 (BRCA RNA), GSE76938 (PRAD methylation), GSE29079 (PRAD RNA).

## Repository contents
- One folder per analysis step (see "Run order"); `archive/` folders hold superseded attempts.
- `docs/script_manifest.md`: provenance and status of every script (how each was reached and whether it was run).
- `CITATION.cff`: citation metadata.
- `LICENSE`: MIT licence for the code and documents in this repository (copyright Abeeb Babatunde Jamiu). It does not cover the TCGA and GEO data, which are not stored here and keep their own terms.
- `.gitignore` excludes raw, processed and result data folders and R session files.

## Things a reader needs to know
- **Keep `GDCdata/`.** Step 7 reads the per-sample RPPA files in `GDCdata/TCGA-BRCA/Proteome_Profiling/Protein_Expression_Quantification` directly (`07_0`, `brca_rppa_dir`). Step 1.4 recreates them.
- **BRCA protein is RPPA only.** 30 patients with paired tumour and normal samples, selected in Step 1.4 as early-stage patients with both sample types.
- **Early-stage cohort definition (Step 5.0).** Stage comes from the RNA metadata, which has no stage for 52 patients, so the reported cohorts (16 three-omics and 51 methylation-plus-RNA pairs) leave out early-stage patients that the clinical file would add (18 and 56). The as-run cohorts are the reported ones. Folder 12 takes stage from the clinical file and uses 56 breast pairs, so its cohorts are defined differently.
- **PRAD manifest.** The Step 0.5 script as run uses AJCC stage for PRAD (eligible 0); the PRAD cohort in Step 5.0 does not use it and defines Early as Gleason score of 7 or below.
- **Re-analyses after Step 10.** Folder 11 addresses three limitations of the original Steps 5, 7 and 8, which the Step 5, 7 and 8 results and the thesis Limitations still describe as they were run: prostate RNA on the count scale in MOFA2 and DIABLO (re-run on the log2 scale), internal validation without repeating feature selection (nested resampling added), and external AUCs without confidence intervals (added).
- **`needs_further_validation` (Step 10 panel tables).** The column is produced by `10_1_final_panel_tables_EDITED.R` and flags only features that could not be evaluated externally (VAV1 and CHD9NB in BRCA). The tables were not regenerated after the nested analysis in folder 11, in which the locked features were not reliably re-selected; the thesis therefore treats all panel features as requiring validation in independent patients.
- **Metabolomics was excluded**; `01_5` records the decision.
- **Vital-status imbalance.** The RPPA patients are 53% dead against 8.5% of the others (Fisher p = 1.6e-09); no survival variable was used in the selection.
- The enrichment queries to KEGG were made online; record the date and package versions with the results.
