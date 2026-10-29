# ============================================================
# STEP 0.5 — Clinical Metadata & Sample Harmonization
# MyThesis_MultiOmics Pipeline
#
# Purpose:
#   Build a single "master sample manifest" that harmonizes
#   TCGA barcodes across epigenomics (1.0), transcriptomics (2.0),
#   and proteomics (3.0) layers, joined to clinical covariates.
#
# Scope boundary (important):
#   - This step is READ-ONLY with respect to 1.0-4.0. It never
#     rewrites, re-normalizes, or re-batch-corrects anything.
#   - It does NOT touch batch/plate correction — that stays inside
#     whatever normalization you already did in 1.0-3.0.
#   - Its only output is a manifest table used later, at 5.1, to
#     FILTER (not re-analyze) your existing result tables.
# ============================================================

## 0. Setup ---------------------------------------------------
library(TCGAbiolinks)
library(dplyr)
library(readr)
library(tibble)

out_dir <- "00.5_clinical_metadata_harmonization"
dir.create(out_dir, showWarnings = FALSE)

## 1. Pull sample IDs from EXISTING outputs (read-only) --------
# Adjust these paths/objects to match your actual 1.0-3.0 outputs.
# The goal is just to get a vector of sample barcodes from each layer.

epi_samples  <- colnames(readRDS("01_epigenomics/dmp_results.rds"))
txn_samples  <- colnames(readRDS("02_transcriptomics/deg_results.rds"))
prot_samples <- colnames(read_csv("03_proteomics/protein_matrix.csv"))

# TCGA barcodes are sample-level (e.g. TCGA-XX-XXXX-01A-...);
# truncate to the 12-character patient-level ID so the three
# layers can be matched to the same person.
to_patient_id <- function(x) substr(x, 1, 12)

epi_ids  <- unique(to_patient_id(epi_samples))
txn_ids  <- unique(to_patient_id(txn_samples))
prot_ids <- unique(to_patient_id(prot_samples))

## 2. Pull clinical data ----------------------------------------
clin_brca <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
clin_prad <- GDCquery_clinic(project = "TCGA-PRAD", type = "clinical")

clinical <- bind_rows(
  clin_brca %>% mutate(project = "BRCA"),
  clin_prad %>% mutate(project = "PRAD")
) %>%
  rename(patient_id = submitter_id) %>%
  select(
    patient_id, project, age_at_index, gender,
    ajcc_pathologic_stage, tumor_grade,
    vital_status, days_to_death, days_to_last_follow_up
  )
  # Add project-specific covariates if present in your clinical pull,
  # e.g. PSA (PRAD) or ER/PR/HER2 status (BRCA).

## 3. Build the crosswalk / manifest -----------------------------
manifest <- tibble(
  patient_id = union(union(epi_ids, txn_ids), prot_ids)
) %>%
  mutate(
    has_epigenomics      = patient_id %in% epi_ids,
    has_transcriptomics  = patient_id %in% txn_ids,
    has_proteomics       = patient_id %in% prot_ids
  ) %>%
  left_join(clinical, by = "patient_id")

## 4. Flag eligibility for MOFA / DIABLO --------------------------
manifest <- manifest %>%
  mutate(
    complete_omics = has_epigenomics & has_transcriptomics & has_proteomics,
    complete_clinical = !is.na(ajcc_pathologic_stage) & !is.na(vital_status),
    eligible_for_integration = complete_omics & complete_clinical
  )

## 5. Attrition report (this is your defensibility trail) ---------
attrition <- tibble(
  stage = c(
    "Total unique patients across any omics layer",
    "Have epigenomics",
    "Have transcriptomics",
    "Have proteomics",
    "Have all 3 omics layers",
    "Have complete clinical annotation",
    "Eligible for MOFA/DIABLO (complete omics + clinical)"
  ),
  n = c(
    nrow(manifest),
    sum(manifest$has_epigenomics),
    sum(manifest$has_transcriptomics),
    sum(manifest$has_proteomics),
    sum(manifest$complete_omics),
    sum(manifest$complete_clinical),
    sum(manifest$eligible_for_integration)
  )
)

## 6. Write outputs (does NOT modify 1.0-3.0 files) ----------------
write_csv(manifest, file.path(out_dir, "sample_manifest_harmonized.csv"))
write_csv(attrition, file.path(out_dir, "attrition_report.csv"))

cat("Step 0.5 complete.\n")
print(attrition)


# ============================================================
# USAGE AT STEP 5.1 (Cross-omics overlap) — subset, don't rerun
# ============================================================
#
# manifest <- read_csv("00.5_clinical_metadata_harmonization/sample_manifest_harmonized.csv")
# eligible_ids <- manifest %>% filter(eligible_for_integration) %>% pull(patient_id)
#
# # Filter your EXISTING result matrices down to the eligible set.
# # This is a column-subset operation, not a re-analysis.
# epi_filtered  <- epi_matrix[,  substr(colnames(epi_matrix), 1, 12)  %in% eligible_ids]
# txn_filtered  <- txn_matrix[,  substr(colnames(txn_matrix), 1, 12)  %in% eligible_ids]
# prot_filtered <- prot_matrix[, substr(colnames(prot_matrix), 1, 12) %in% eligible_ids]
#
# # These filtered objects — not the originals — feed into 5.2 (MOFA)
# # and 5.3 (DIABLO). Your 1.0-4.0 result tables remain untouched.