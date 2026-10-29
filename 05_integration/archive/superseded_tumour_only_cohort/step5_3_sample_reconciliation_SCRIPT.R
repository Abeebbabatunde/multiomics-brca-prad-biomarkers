# ============================================================
# STEP 5.3 — FINAL SAMPLE RECONCILIATION
# ============================================================
#
# Purpose:
#   1. Identify BRCA patients without primary-tumor methylation
#   2. Identify PRAD patients with multiple primary-tumor samples
#   3. Inspect complete TCGA barcodes, including vial information
#   4. Reconcile BRCA protein identifiers
#   5. Create preliminary final sample manifests
#
# IMPORTANT:
#   NO original RNA, methylation, or protein matrix is altered.
#   This step only creates diagnostic tables/manifests.
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP 5.3 — FINAL SAMPLE RECONCILIATION\n")
cat("============================================================\n")


# ============================================================
# 1. OUTPUT DIRECTORY
# ============================================================

dir.create(
  "05_FinalResults/Step5_SampleMaps",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. HELPER FUNCTION — FULL TCGA SAMPLE INFORMATION
# ============================================================

tcga_sample_code <- function(x) {

  x <- as.character(x)

  out <- sub(
    "^TCGA-[^-]+-[^-]+-([0-9]{2}).*",
    "\\1",
    x
  )

  out[
    !grepl(
      "^TCGA-[^-]+-[^-]+-[0-9]{2}",
      x
    )
  ] <- NA

  out
}


tcga_vial <- function(x) {

  x <- as.character(x)

  out <- sub(
    "^TCGA-[^-]+-[^-]+-[0-9]{2}([A-Z]).*",
    "\\1",
    x
  )

  out[
    !grepl(
      "^TCGA-[^-]+-[^-]+-[0-9]{2}[A-Z]",
      x
    )
  ] <- NA

  out
}


# ============================================================
# 3. ADD FULL SAMPLE INFORMATION TO EXISTING MAPS
# ============================================================

brca_rna_map$sample_code_full <- tcga_sample_code(
  brca_rna_map$sample_id
)

brca_rna_map$vial <- tcga_vial(
  brca_rna_map$sample_id
)


brca_meth_map$sample_code_full <- tcga_sample_code(
  brca_meth_map$sample_id
)

brca_meth_map$vial <- tcga_vial(
  brca_meth_map$sample_id
)


prad_rna_map$sample_code_full <- tcga_sample_code(
  prad_rna_map$sample_id
)

prad_rna_map$vial <- tcga_vial(
  prad_rna_map$sample_id
)


prad_meth_map$sample_code_full <- tcga_sample_code(
  prad_meth_map$sample_id
)

prad_meth_map$vial <- tcga_vial(
  prad_meth_map$sample_id
)


# ============================================================
# PART A — BRCA PRIMARY-TUMOR METHYLATION AUDIT
# ============================================================

cat("\n============================================================\n")
cat("5.3A — BRCA PRIMARY-TUMOR METHYLATION AUDIT\n")
cat("============================================================\n")


brca_meth_01 <- brca_meth_map[
  brca_meth_map$target_patient &
    brca_meth_map$sample_code_full == "01",
]


brca_meth_01_patients <- unique(
  brca_meth_01$patient_id
)


brca_missing_meth_01 <- setdiff(
  brca_ids,
  brca_meth_01_patients
)


cat(
  "\nBRCA target patients:",
  length(brca_ids),
  "\n"
)

cat(
  "BRCA patients with methylation 01:",
  length(brca_meth_01_patients),
  "\n"
)

cat(
  "BRCA patients without methylation 01:",
  length(brca_missing_meth_01),
  "\n"
)


cat("\nPatients without methylation 01:\n")

if (length(brca_missing_meth_01) > 0) {

  print(brca_missing_meth_01)

} else {

  cat("NONE\n")

}


# Save missing-primary-tumor table

brca_missing_meth_table <- data.frame(
  patient_id = brca_missing_meth_01,
  stringsAsFactors = FALSE
)


write.csv(
  brca_missing_meth_table,
  "05_FinalResults/Step5_SampleMaps/BRCA_missing_primary_methylation.csv",
  row.names = FALSE
)


# ============================================================
# 4. SHOW OTHER BRCA METHYLATION SAMPLES FOR MISSING PATIENTS
# ============================================================

cat("\nBRCA methylation samples for patients lacking 01:\n")


brca_missing_meth_details <- brca_meth_map[
  brca_meth_map$patient_id %in% brca_missing_meth_01,
]


if (nrow(brca_missing_meth_details) > 0) {

  print(
    brca_missing_meth_details[
      ,
      c(
        "patient_id",
        "sample_id",
        "sample_code_full",
        "vial"
      )
    ]
  )

} else {

  cat("No methylation samples found for these patients.\n")

}


write.csv(
  brca_missing_meth_details,
  "05_FinalResults/Step5_SampleMaps/BRCA_missing_primary_methylation_details.csv",
  row.names = FALSE
)


# ============================================================
# PART B — PRAD DUPLICATE PRIMARY-TUMOR AUDIT
# ============================================================

cat("\n============================================================\n")
cat("5.3B — PRAD DUPLICATE PRIMARY-TUMOR AUDIT\n")
cat("============================================================\n")


# ------------------------------------------------------------
# PRAD RNA
# ------------------------------------------------------------

prad_rna_01 <- prad_rna_map[
  prad_rna_map$target_patient &
    prad_rna_map$sample_code_full == "01",
]


prad_rna_01_counts <- table(
  prad_rna_01$patient_id
)


prad_rna_duplicate_patients <- names(
  prad_rna_01_counts[
    prad_rna_01_counts > 1
  ]
)


cat(
  "\nPRAD RNA patients with >1 primary tumor sample:",
  length(prad_rna_duplicate_patients),
  "\n"
)


if (length(prad_rna_duplicate_patients) > 0) {

  print(prad_rna_duplicate_patients)

}


prad_rna_duplicate_details <- prad_rna_01[
  prad_rna_01$patient_id %in%
    prad_rna_duplicate_patients,
]


cat("\nPRAD RNA duplicate primary samples:\n")

if (nrow(prad_rna_duplicate_details) > 0) {

  print(
    prad_rna_duplicate_details[
      ,
      c(
        "patient_id",
        "sample_id",
        "sample_code_full",
        "vial"
      )
    ]
  )

} else {

  cat("NONE\n")

}


# ------------------------------------------------------------
# PRAD METHYLATION
# ------------------------------------------------------------

prad_meth_01 <- prad_meth_map[
  prad_meth_map$target_patient &
    prad_meth_map$sample_code_full == "01",
]


prad_meth_01_counts <- table(
  prad_meth_01$patient_id
)


prad_meth_duplicate_patients <- names(
  prad_meth_01_counts[
    prad_meth_01_counts > 1
  ]
)


cat(
  "\nPRAD methylation patients with >1 primary tumor sample:",
  length(prad_meth_duplicate_patients),
  "\n"
)


if (length(prad_meth_duplicate_patients) > 0) {

  print(prad_meth_duplicate_patients)

}


prad_meth_duplicate_details <- prad_meth_01[
  prad_meth_01$patient_id %in%
    prad_meth_duplicate_patients,
]


cat("\nPRAD methylation duplicate primary samples:\n")

if (nrow(prad_meth_duplicate_details) > 0) {

  print(
    prad_meth_duplicate_details[
      ,
      c(
        "patient_id",
        "sample_id",
        "sample_code_full",
        "vial"
      )
    ]
  )

} else {

  cat("NONE\n")

}


# ============================================================
# 5. COMBINE PRAD DUPLICATE INFORMATION
# ============================================================

prad_duplicate_details <- rbind(

  data.frame(
    omics = "RNA",
    prad_rna_duplicate_details,
    stringsAsFactors = FALSE
  ),

  data.frame(
    omics = "Methylation",
    prad_meth_duplicate_details,
    stringsAsFactors = FALSE
  )

)


write.csv(
  prad_duplicate_details,
  "05_FinalResults/Step5_SampleMaps/PRAD_duplicate_primary_samples.csv",
  row.names = FALSE
)


# ============================================================
# PART C — BRCA RNA/METHYLATION PRIMARY-TUMOR MATCHING
# ============================================================

cat("\n============================================================\n")
cat("5.3C — BRCA PRIMARY-TUMOR RNA/METHYLATION MATCHING\n")
cat("============================================================\n")


brca_rna_01_patients <- unique(
  brca_rna_map$patient_id[
    brca_rna_map$target_patient &
      brca_rna_map$sample_code_full == "01"
  ]
)


brca_meth_01_patients <- unique(
  brca_meth_map$patient_id[
    brca_meth_map$target_patient &
      brca_meth_map$sample_code_full == "01"
  ]
)


brca_rna_meth_01_common <- intersect(
  brca_rna_01_patients,
  brca_meth_01_patients
)


brca_rna_only_01 <- setdiff(
  brca_rna_01_patients,
  brca_meth_01_patients
)


brca_meth_only_01 <- setdiff(
  brca_meth_01_patients,
  brca_rna_01_patients
)


cat(
  "\nBRCA RNA primary-tumor patients:",
  length(brca_rna_01_patients),
  "\n"
)

cat(
  "BRCA methylation primary-tumor patients:",
  length(brca_meth_01_patients),
  "\n"
)

cat(
  "BRCA matched primary-tumor patients:",
  length(brca_rna_meth_01_common),
  "\n"
)

cat(
  "BRCA RNA-only primary-tumor patients:",
  length(brca_rna_only_01),
  "\n"
)

cat(
  "BRCA methylation-only primary-tumor patients:",
  length(brca_meth_only_01),
  "\n"
)


# ============================================================
# PART D — PRAD RNA/METHYLATION PRIMARY-TUMOR MATCHING
# ============================================================

cat("\n============================================================\n")
cat("5.3D — PRAD PRIMARY-TUMOR RNA/METHYLATION MATCHING\n")
cat("============================================================\n")


prad_rna_01_patients <- unique(
  prad_rna_map$patient_id[
    prad_rna_map$target_patient &
      prad_rna_map$sample_code_full == "01"
  ]
)


prad_meth_01_patients <- unique(
  prad_meth_map$patient_id[
    prad_meth_map$target_patient &
      prad_meth_map$sample_code_full == "01"
  ]
)


prad_rna_meth_01_common <- intersect(
  prad_rna_01_patients,
  prad_meth_01_patients
)


cat(
  "\nPRAD RNA primary-tumor patients:",
  length(prad_rna_01_patients),
  "\n"
)

cat(
  "PRAD methylation primary-tumor patients:",
  length(prad_meth_01_patients),
  "\n"
)

cat(
  "PRAD matched primary-tumor patients:",
  length(prad_rna_meth_01_common),
  "\n"
)


# ============================================================
# PART E — PRELIMINARY BRCA FINAL SAMPLE MANIFEST
# ============================================================
#
# IMPORTANT:
#   We do NOT yet select between duplicate samples.
#   This manifest records all available target samples.
# ============================================================

cat("\n============================================================\n")
cat("5.3E — BRCA PRELIMINARY SAMPLE MANIFEST\n")
cat("============================================================\n")


brca_manifest <- data.frame(
  patient_id = brca_ids,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# RNA primary tumor
# ------------------------------------------------------------

brca_rna_01_lookup <- brca_rna_map[
  brca_rna_map$target_patient &
    brca_rna_map$sample_code_full == "01",
  c(
    "patient_id",
    "sample_id",
    "vial"
  )
]


colnames(brca_rna_01_lookup) <- c(
  "patient_id",
  "RNA_01_sample",
  "RNA_01_vial"
)


brca_manifest <- merge(
  brca_manifest,
  brca_rna_01_lookup,
  by = "patient_id",
  all.x = TRUE,
  sort = FALSE
)


# ------------------------------------------------------------
# Methylation primary tumor
# ------------------------------------------------------------

brca_meth_01_lookup <- brca_meth_map[
  brca_meth_map$target_patient &
    brca_meth_map$sample_code_full == "01",
  c(
    "patient_id",
    "sample_id",
    "vial"
  )
]


colnames(brca_meth_01_lookup) <- c(
  "patient_id",
  "METH_01_sample",
  "METH_01_vial"
)


brca_manifest <- merge(
  brca_manifest,
  brca_meth_01_lookup,
  by = "patient_id",
  all.x = TRUE,
  sort = FALSE
)


# Restore target patient order

brca_manifest <- brca_manifest[
  match(
    brca_ids,
    brca_manifest$patient_id
  ),
]


# ------------------------------------------------------------
# Match status
# ------------------------------------------------------------

brca_manifest$RNA_METH_01_matched <-

  !is.na(brca_manifest$RNA_01_sample) &
  !is.na(brca_manifest$METH_01_sample)


print(head(brca_manifest, 20))


write.csv(
  brca_manifest,
  "05_FinalResults/Step5_SampleMaps/BRCA_preliminary_sample_manifest.csv",
  row.names = FALSE
)


# ============================================================
# PART F — PRELIMINARY PRAD FINAL SAMPLE MANIFEST
# ============================================================

cat("\n============================================================\n")
cat("5.3F — PRAD PRELIMINARY SAMPLE MANIFEST\n")
cat("============================================================\n")


prad_manifest <- data.frame(
  patient_id = prad_ids,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# RNA primary tumor
# ------------------------------------------------------------

prad_rna_01_lookup <- prad_rna_map[
  prad_rna_map$target_patient &
    prad_rna_map$sample_code_full == "01",
  c(
    "patient_id",
    "sample_id",
    "vial"
  )
]


colnames(prad_rna_01_lookup) <- c(
  "patient_id",
  "RNA_01_sample",
  "RNA_01_vial"
)


prad_manifest <- merge(
  prad_manifest,
  prad_rna_01_lookup,
  by = "patient_id",
  all.x = TRUE,
  sort = FALSE
)


# ------------------------------------------------------------
# Methylation primary tumor
# ------------------------------------------------------------

prad_meth_01_lookup <- prad_meth_map[
  prad_meth_map$target_patient &
    prad_meth_map$sample_code_full == "01",
  c(
    "patient_id",
    "sample_id",
    "vial"
  )
]


colnames(prad_meth_01_lookup) <- c(
  "patient_id",
  "METH_01_sample",
  "METH_01_vial"
)


prad_manifest <- merge(
  prad_manifest,
  prad_meth_01_lookup,
  by = "patient_id",
  all.x = TRUE,
  sort = FALSE
)


# Restore target patient order

prad_manifest <- prad_manifest[
  match(
    prad_ids,
    prad_manifest$patient_id
  ),
]


# ------------------------------------------------------------
# Match status
# ------------------------------------------------------------

prad_manifest$RNA_METH_01_matched <-

  !is.na(prad_manifest$RNA_01_sample) &
  !is.na(prad_manifest$METH_01_sample)


print(head(prad_manifest, 20))


write.csv(
  prad_manifest,
  "05_FinalResults/Step5_SampleMaps/PRAD_preliminary_sample_manifest.csv",
  row.names = FALSE
)


# ============================================================
# PART G — FINAL SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP 5.3 DIAGNOSTIC RECONCILIATION COMPLETE\n")
cat("============================================================\n")


cat("\n---------------- BRCA ----------------\n")

cat(
  "Target patients:",
  length(brca_ids),
  "\n"
)

cat(
  "RNA 01 patients:",
  length(brca_rna_01_patients),
  "\n"
)

cat(
  "Methylation 01 patients:",
  length(brca_meth_01_patients),
  "\n"
)

cat(
  "Matched RNA + methylation 01 patients:",
  length(brca_rna_meth_01_common),
  "\n"
)

cat(
  "Missing methylation 01 patients:",
  length(brca_missing_meth_01),
  "\n"
)


cat("\n---------------- PRAD ----------------\n")

cat(
  "Target patients:",
  length(prad_ids),
  "\n"
)

cat(
  "RNA 01 patients:",
  length(prad_rna_01_patients),
  "\n"
)

cat(
  "Methylation 01 patients:",
  length(prad_meth_01_patients),
  "\n"
)

cat(
  "Matched RNA + methylation 01 patients:",
  length(prad_rna_meth_01_common),
  "\n"
)

cat(
  "PRAD RNA duplicate patients:",
  length(prad_rna_duplicate_patients),
  "\n"
)

cat(
  "PRAD methylation duplicate patients:",
  length(prad_meth_duplicate_patients),
  "\n"
)


cat("\n------------------------------------------------------------\n")
cat("IMPORTANT: No original matrices were modified.\n")
cat("No samples were deleted.\n")
cat("No duplicate samples were selected or discarded.\n")
cat("------------------------------------------------------------\n")