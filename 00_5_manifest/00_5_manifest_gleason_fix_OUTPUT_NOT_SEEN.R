# ============================================================
# STEP 0.5 — Clinical Metadata & Sample Harmonization
# MyThesis_MultiOmics Pipeline
#
# Design decision (confirmed):
#   - BRCA integration uses 3 omics blocks: epigenomics, transcriptomics,
#     proteomics (CPTAC/TCGA-linked, barcodes confirmed compatible).
#   - PRAD integration uses 2 omics blocks: epigenomics, transcriptomics.
#     PRAD proteomics (E-PROT-21 / Iglesias-Gato cohort) is EXCLUDED from
#     the manifest and from MOFA/DIABLO — it is an independent cohort
#     with no shared patient identity to TCGA-PRAD, confirmed via its
#     experiment-design file (no TCGA barcode, no TCGA reference).
#     E-PROT-21 is reserved for a separate downstream validation step
#     (checking direction of change for candidate biomarkers), not for
#     integration.
# ============================================================

library(dplyr)
library(readr)
library(tibble)

out_dir <- "00.5_clinical_metadata_harmonization"
dir.create(out_dir, showWarnings = FALSE)

to_patient_id <- function(x) substr(x, 1, 12)

## ------------------------------------------------------------
## PART A — Load sample-level matrices
## ------------------------------------------------------------

meth_brca <- read_csv("02_ProcessedData/Breast/brca_meth_normalized.csv", show_col_types = FALSE)
rna_brca  <- readRDS("02_ProcessedData/Breast/dds_brca.rds")
prot_brca <- read_csv("02_ProcessedData/Breast/brca_protein_normalized.csv", show_col_types = FALSE)

meth_prad <- read_csv("02_ProcessedData/Prostate/prad_meth_normalized.csv", show_col_types = FALSE)
rna_prad  <- readRDS("02_ProcessedData/Prostate/dds_prad.rds")
# PRAD proteomics intentionally NOT loaded here — see design note above.

## ------------------------------------------------------------
## PART B — Extract patient-level IDs (first 12 chars of TCGA barcode)
## ------------------------------------------------------------

epi_ids_brca  <- unique(to_patient_id(colnames(meth_brca)[-1]))          # drop "...1" feature-id col
txn_ids_brca  <- unique(to_patient_id(colnames(rna_brca)))
prot_ids_brca <- unique(to_patient_id(gsub("\\.", "-", colnames(prot_brca)[-1])))  # dots -> dashes

epi_ids_prad  <- unique(to_patient_id(colnames(meth_prad)[-1]))
txn_ids_prad  <- unique(to_patient_id(colnames(rna_prad)))

## ------------------------------------------------------------
## PART C — Clinical data (from your existing local files)
## ------------------------------------------------------------

clin_brca <- read_csv("01_RawData/Breast/clinical_brca.csv", show_col_types = FALSE) %>%
  transmute(
    patient_id = to_patient_id(submitter_id),
    project = "BRCA",
    ajcc_pathologic_stage, tumor_grade, vital_status,
    days_to_death, days_to_last_follow_up
  )

clin_prad <- read_csv("01_RawData/Prostate/clinical_prad.csv", show_col_types = FALSE) %>%
  transmute(
    patient_id = to_patient_id(submitter_id),
    project = "PRAD",
    ajcc_pathologic_stage, tumor_grade, vital_status,
    gleason_score,
    days_to_death, days_to_last_follow_up
  )

clinical <- bind_rows(clin_brca, clin_prad) %>% distinct(patient_id, .keep_all = TRUE)

## ------------------------------------------------------------
## PART D — Build manifests
## BRCA: 3-omics eligibility gate.
## PRAD: 2-omics eligibility gate (no proteomics column at all).
## ------------------------------------------------------------

manifest_brca <- tibble(patient_id = union(union(epi_ids_brca, txn_ids_brca), prot_ids_brca)) %>%
  mutate(
    has_epigenomics     = patient_id %in% epi_ids_brca,
    has_transcriptomics = patient_id %in% txn_ids_brca,
    has_proteomics      = patient_id %in% prot_ids_brca
  ) %>%
  left_join(clinical %>% filter(project == "BRCA"), by = "patient_id") %>%
  mutate(
    complete_omics    = has_epigenomics & has_transcriptomics & has_proteomics,
    complete_clinical = !is.na(vital_status) & !is.na(ajcc_pathologic_stage),
    eligible_for_integration = complete_omics & complete_clinical
  )

manifest_prad <- tibble(patient_id = union(epi_ids_prad, txn_ids_prad)) %>%
  mutate(
    has_epigenomics     = patient_id %in% epi_ids_prad,
    has_transcriptomics = patient_id %in% txn_ids_prad
    # no has_proteomics column for PRAD — proteomics is not part of
    # this cohort's integration (see design note above)
  ) %>%
  left_join(clinical %>% filter(project == "PRAD"), by = "patient_id") %>%
  mutate(
    complete_omics    = has_epigenomics & has_transcriptomics,
    # PRAD is staged by Gleason score, not AJCC pathologic stage
    # (that field is largely unpopulated for TCGA-PRAD) — use
    # gleason_score as the staging criterion here instead.
    complete_clinical = !is.na(vital_status) & !is.na(gleason_score),
    eligible_for_integration = complete_omics & complete_clinical
  )

manifest <- bind_rows(manifest_brca, manifest_prad)

## ------------------------------------------------------------
## PART E — Attrition report + write outputs
## ------------------------------------------------------------

attrition <- manifest %>%
  group_by(project) %>%
  summarise(
    total = n(),
    has_epi = sum(has_epigenomics, na.rm = TRUE),
    has_txn = sum(has_transcriptomics, na.rm = TRUE),
    has_prot = ifelse("has_proteomics" %in% names(manifest), sum(has_proteomics, na.rm = TRUE), NA),
    complete_omics = sum(complete_omics, na.rm = TRUE),
    complete_clinical = sum(complete_clinical, na.rm = TRUE),
    eligible_for_integration = sum(eligible_for_integration, na.rm = TRUE)
  )

write_csv(manifest, file.path(out_dir, "sample_manifest_harmonized.csv"))
write_csv(attrition, file.path(out_dir, "attrition_report.csv"))

cat("Step 0.5 complete.\n")
cat("BRCA: 3-omics integration (epi + txn + prot)\n")
cat("PRAD: 2-omics integration (epi + txn); E-PROT-21 reserved for validation only\n\n")
print(attrition)


# ============================================================
# USAGE AT STEP 5.1 (Cross-omics overlap) — subset, don't rerun
# ============================================================
#
# manifest <- read_csv("00.5_clinical_metadata_harmonization/sample_manifest_harmonized.csv")
#
# brca_eligible <- manifest %>% filter(project == "BRCA", eligible_for_integration) %>% pull(patient_id)
# prad_eligible <- manifest %>% filter(project == "PRAD", eligible_for_integration) %>% pull(patient_id)
#
# # BRCA: subset all 3 blocks
# epi_brca_filtered  <- meth_brca[,  c(1, which(substr(colnames(meth_brca)[-1], 1, 12)  %in% brca_eligible) + 1)]
# txn_brca_filtered  <- rna_brca[,  substr(colnames(rna_brca), 1, 12)  %in% brca_eligible]
# prot_brca_filtered <- prot_brca[, c(1, which(substr(gsub("\\.", "-", colnames(prot_brca)[-1]), 1, 12) %in% brca_eligible) + 1)]
#
# # PRAD: subset only 2 blocks
# epi_prad_filtered  <- meth_prad[, c(1, which(substr(colnames(meth_prad)[-1], 1, 12) %in% prad_eligible) + 1)]
# txn_prad_filtered  <- rna_prad[,  substr(colnames(rna_prad), 1, 12) %in% prad_eligible]
#
# # PRAD proteomics (E-PROT-21) is NOT subset here — it feeds a separate
# # validation step after biomarker candidates are identified, not 5.1-5.3.