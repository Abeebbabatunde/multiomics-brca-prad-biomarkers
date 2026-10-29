# MyThesis_MultiOmics: script library

In silico multi-omics biomarker discovery in early-stage breast (TCGA-BRCA) and prostate (TCGA-PRAD) cancer, written in R. Every script here is the code that produced the reported results, filed with a label that says how it reached this library (see "Script labels"). Nothing was written afresh for the library.

## Citation
If you use these scripts, please cite:

Jamiu AB (2026). Multi-omics biomarker discovery in early-stage breast and prostate cancer: analysis scripts. GitHub. https://github.com/Abeebbabatunde/multiomics-brca-prad-biomarkers

`CITATION.cff` holds the same information in machine-readable form (GitHub shows it as "Cite this repository").

## Setup
- R 4.6.0 (2026-04-24), Bioconductor 3.23, Windows 11 (x64).
- Packages: TCGAbiolinks 2.40.0, SummarizedExperiment 1.42.0, minfi 1.58.0, limma 3.68.2, edgeR 4.10.0, DESeq2 1.52.0, MOFA2 1.22.0, mixOmics 6.36.0, glmnet 5.0, randomForest 4.7.1.2, e1071 1.7.17, survival 3.8.6, sesameData 1.30.0. Step 7 also needs kernlab (SVM-RFE). The pre-submission re-analyses (`11_pre_submission_reanalyses/`) also use Boruta, pROC and mofapy2 0.7.5 (Python 3.11); their versions (Boruta 10.0.0, pROC 1.19.0.1) and the full `sessionInfo()` are in `11_pre_submission_reanalyses/`.
- Seeds: 123 (MOFA2, Step 5), 1 (shuffled-label check, Step 5), 42 (Steps 7, 8 and 9.1), 2026 (all re-analyses in folder 11). The Step 4 GSEA scripts set no seed.
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

Each step folder lists its scripts in numeric order (`05_0a`, `05_0b`, `05_1a` and so on).

## Script labels (end of file name)
- `AS_RUN`: the script as it was run, matched to the printed results.
- `AS_PASTED`: pasted from the session and filed unchanged.
- `EDITED`: changed after a first run; only the protein source lines differ (BRCA protein is RPPA throughout). The first-run version is in that folder's `archive/`.
- `RECOVERED_FROM_HISTORY`: copied from the RStudio history (20 Sep 2026).
- `NOT_CONFIRMED_RUN`, `OUTPUT_NOT_SEEN`: present for completeness; their output was not confirmed.
- No label: reconciled against the printed results and the document.
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
- **Early-stage cohort definition (Step 5.0).** Stage comes from the RNA metadata, which has no stage for 52 patients, so the reported cohorts (16 three-omics and 51 methylation-plus-RNA pairs) leave out early-stage patients that the clinical file would add (18 and 56). The as-run cohorts are the reported ones.
- **PRAD manifest.** The Step 0.5 script as run uses AJCC stage for PRAD (eligible 0); the PRAD cohort in Step 5.0 does not use it and defines Early as Gleason score of 7 or below.
- **Re-analyses after Step 10.** Folder 11 addresses three limitations of the original Steps 5, 7 and 8, which the Step 5, 7 and 8 results and the thesis Limitations still describe as they were run: prostate RNA on the count scale in MOFA2 and DIABLO (re-run on the log2 scale), internal validation without repeating feature selection (nested resampling added), and external AUCs without confidence intervals (added).
- **`needs_further_validation` (Step 10 panel tables).** The column is produced by `10_1_final_panel_tables_EDITED.R` and flags only features that could not be evaluated externally (VAV1 and CHD9NB in BRCA). The tables were not regenerated after the nested analysis in folder 11, in which the locked features were not reliably re-selected; the thesis therefore treats all panel features as requiring validation in independent patients.
- **Metabolomics was excluded**; `01_5` records the decision.
- **Vital-status imbalance.** The RPPA patients are 53% dead against 8.5% of the others (Fisher p = 1.6e-09); no survival variable was used in the selection.
- The enrichment queries to KEGG were made online; record the date and package versions with the results.
