# ============================================================
# STEP 5.1 — BUILD FINAL MULTI-OMICS INTEGRATION COHORTS
# ============================================================

cat("\n============================================================\n")
cat("STEP 5.1 — FINAL MULTI-OMICS COHORT CONSTRUCTION\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Load harmonized manifest
# ------------------------------------------------------------

manifest <- read.csv(
  "00.5_clinical_metadata_harmonization/sample_manifest_harmonized.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# PART A — BRCA 3-OMICS COHORT
# ============================================================

cat("\n============================================================\n")
cat("PART A — BRCA 3-OMICS COHORT\n")
cat("============================================================\n")


brca_3omics <- subset(
  manifest,
  project == "BRCA" &
    has_epigenomics == TRUE &
    has_transcriptomics == TRUE &
    has_proteomics == TRUE
)


# ------------------------------------------------------------
# Check BRCA cohort
# ------------------------------------------------------------

cat("\nNumber of BRCA 3-omics patients:",
    nrow(brca_3omics),
    "\n")

cat("\nBRCA 3-omics patient IDs:\n")
print(brca_3omics$patient_id)


# ------------------------------------------------------------
# Check duplicates
# ------------------------------------------------------------

cat("\nDuplicate BRCA patient IDs:\n")

brca_duplicates <- brca_3omics$patient_id[
  duplicated(brca_3omics$patient_id)
]

print(brca_duplicates)

cat(
  "Number of duplicated BRCA IDs:",
  length(brca_duplicates),
  "\n"
)


# ------------------------------------------------------------
# Verify all three flags
# ------------------------------------------------------------

cat("\nBRCA omics flags:\n")

print(
  table(
    brca_3omics$has_epigenomics,
    brca_3omics$has_transcriptomics,
    brca_3omics$has_proteomics
  )
)


# ------------------------------------------------------------
# Save BRCA cohort
# ------------------------------------------------------------

write.csv(
  brca_3omics,
  "05_FinalResults/Step5_BRCA_3omics_cohort.csv",
  row.names = FALSE
)

saveRDS(
  brca_3omics,
  "05_FinalResults/Step5_BRCA_3omics_cohort.rds"
)


# ============================================================
# PART B — PRAD 2-OMICS COHORT
# ============================================================

cat("\n============================================================\n")
cat("PART B — PRAD 2-OMICS COHORT\n")
cat("============================================================\n")


prad_2omics <- subset(
  manifest,
  project == "PRAD" &
    has_epigenomics == TRUE &
    has_transcriptomics == TRUE
)


# ------------------------------------------------------------
# Check PRAD cohort
# ------------------------------------------------------------

cat("\nNumber of PRAD 2-omics patients:",
    nrow(prad_2omics),
    "\n")

cat("\nFirst PRAD 2-omics patient IDs:\n")
print(head(prad_2omics$patient_id, 20))


# ------------------------------------------------------------
# Check duplicates
# ------------------------------------------------------------

cat("\nDuplicate PRAD patient IDs:\n")

prad_duplicates <- prad_2omics$patient_id[
  duplicated(prad_2omics$patient_id)
]

print(prad_duplicates)

cat(
  "Number of duplicated PRAD IDs:",
  length(prad_duplicates),
  "\n"
)


# ------------------------------------------------------------
# Verify two-omics flags
# ------------------------------------------------------------

cat("\nPRAD omics flags:\n")

print(
  table(
    prad_2omics$has_epigenomics,
    prad_2omics$has_transcriptomics
  )
)


# ------------------------------------------------------------
# Save PRAD cohort
# ------------------------------------------------------------

write.csv(
  prad_2omics,
  "05_FinalResults/Step5_PRAD_2omics_cohort.csv",
  row.names = FALSE
)

saveRDS(
  prad_2omics,
  "05_FinalResults/Step5_PRAD_2omics_cohort.rds"
)


# ============================================================
# PART C — FINAL COHORT SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("FINAL COHORT SUMMARY\n")
cat("============================================================\n")


final_cohort_summary <- data.frame(
  Cancer = c(
    "BRCA",
    "PRAD"
  ),

  Integration_Type = c(
    "3-omics",
    "2-omics"
  ),

  Epigenomics = c(
    TRUE,
    TRUE
  ),

  Transcriptomics = c(
    TRUE,
    TRUE
  ),

  Proteomics = c(
    TRUE,
    FALSE
  ),

  N_patients = c(
    nrow(brca_3omics),
    nrow(prad_2omics)
  ),

  stringsAsFactors = FALSE
)


print(final_cohort_summary)


# ------------------------------------------------------------
# Save summary
# ------------------------------------------------------------

write.csv(
  final_cohort_summary,
  "05_FinalResults/Step5_final_integration_cohort_summary.csv",
  row.names = FALSE
)

saveRDS(
  final_cohort_summary,
  "05_FinalResults/Step5_final_integration_cohort_summary.rds"
)


# ============================================================
# PART D — SAVE JUST THE PATIENT IDs
# ============================================================

brca_3omics_ids <- brca_3omics$patient_id

prad_2omics_ids <- prad_2omics$patient_id


saveRDS(
  brca_3omics_ids,
  "05_FinalResults/Step5_BRCA_3omics_patient_ids.rds"
)

saveRDS(
  prad_2omics_ids,
  "05_FinalResults/Step5_PRAD_2omics_patient_ids.rds"
)


write.csv(
  data.frame(patient_id = brca_3omics_ids),
  "05_FinalResults/Step5_BRCA_3omics_patient_ids.csv",
  row.names = FALSE
)

write.csv(
  data.frame(patient_id = prad_2omics_ids),
  "05_FinalResults/Step5_PRAD_2omics_patient_ids.csv",
  row.names = FALSE
)


# ============================================================
# FINAL CHECK
# ============================================================

cat("\n============================================================\n")
cat("FINAL CHECK\n")
cat("============================================================\n")

cat(
  "\nBRCA 3-omics patients:",
  length(brca_3omics_ids),
  "\n"
)

cat(
  "PRAD 2-omics patients:",
  length(prad_2omics_ids),
  "\n"
)

cat(
  "\nBRCA IDs unique:",
  length(unique(brca_3omics_ids)) == length(brca_3omics_ids),
  "\n"
)

cat(
  "PRAD IDs unique:",
  length(unique(prad_2omics_ids)) == length(prad_2omics_ids),
  "\n"
)


cat("\n============================================================\n")
cat("STEP 5.1 COMPLETE\n")
cat("============================================================\n")