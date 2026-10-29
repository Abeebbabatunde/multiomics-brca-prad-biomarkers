# ============================================================
# STEP 5.0G — AUDIT HARMONIZED SAMPLE IDs
# ============================================================

# ------------------------------------------------------------
# 1. Load harmonized manifest
# ------------------------------------------------------------

manifest <- read.csv(
  "00.5_clinical_metadata_harmonization/sample_manifest_harmonized.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------
# 2. BRCA integration candidates
# ------------------------------------------------------------

brca_manifest <- subset(
  manifest,
  project == "BRCA"
)

cat("\n================ BRCA MANIFEST ================\n")

cat("Total BRCA patients:", nrow(brca_manifest), "\n")

cat("All 3 omics:\n")
print(table(brca_manifest$complete_omics))

cat("Eligible for integration:\n")
print(table(brca_manifest$eligible_for_integration))

cat("\nPatients with all 3 omics:\n")

brca_complete <- brca_manifest[
  brca_manifest$complete_omics == TRUE,
]

print(
  brca_complete[, c(
    "patient_id",
    "has_epigenomics",
    "has_transcriptomics",
    "has_proteomics",
    "complete_omics",
    "eligible_for_integration"
  )]
)


# ------------------------------------------------------------
# 3. PRAD integration candidates
# ------------------------------------------------------------

prad_manifest <- subset(
  manifest,
  project == "PRAD"
)

cat("\n================ PRAD MANIFEST ================\n")

cat("Total PRAD patients:", nrow(prad_manifest), "\n")

cat("Epigenomics:\n")
print(table(prad_manifest$has_epigenomics))

cat("Transcriptomics:\n")
print(table(prad_manifest$has_transcriptomics))

cat("Proteomics:\n")
print(table(prad_manifest$has_proteomics))

cat("Complete omics:\n")
print(table(prad_manifest$complete_omics))

cat("Eligible:\n")
print(table(prad_manifest$eligible_for_integration))


# ------------------------------------------------------------
# 4. Inspect PRAD protein metadata
# ------------------------------------------------------------

cat("\n================ PRAD PROTEIN METADATA ================\n")

prad_prot_meta <- readRDS(
  "03_Results/Proteomics/prad_protein_metadata.rds"
)

print(class(prad_prot_meta))
print(dim(prad_prot_meta))
print(colnames(prad_prot_meta))
print(prad_prot_meta)


# ------------------------------------------------------------
# 5. Inspect PRAD protein matrix sample IDs
# ------------------------------------------------------------

cat("\n================ PRAD PROTEIN MATRIX ================\n")

prad_prot <- readRDS(
  "03_Results/Proteomics/prad_protein_normalized.rds"
)

cat("Dimensions:\n")
print(dim(prad_prot))

cat("\nProtein matrix sample IDs:\n")
print(colnames(prad_prot))


# ------------------------------------------------------------
# 6. Inspect PRAD raw proteomics group information
# ------------------------------------------------------------

cat("\n================ PRAD PROTEOMICS GROUPS ================\n")

prad_prot_groups <- readRDS(
  "01_RawData/Prostate/prad_prot_groups.rds"
)

print(class(prad_prot_groups))
print(length(prad_prot_groups))
print(prad_prot_groups)


# ------------------------------------------------------------
# 7. Compare PRAD protein IDs against PRAD clinical IDs
# ------------------------------------------------------------

cat("\n================ PRAD ID COMPARISON ================\n")

clinical_prad <- read.csv(
  "01_RawData/Prostate/clinical_prad.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Clinical PRAD rows:", nrow(clinical_prad), "\n")

cat("\nClinical ID columns:\n")
print(colnames(clinical_prad)[
  grepl(
    "submitter|patient|barcode|case|sample",
    colnames(clinical_prad),
    ignore.case = TRUE
  )
])

cat("\nFirst clinical IDs:\n")

if ("submitter_id" %in% colnames(clinical_prad)) {
  print(head(clinical_prad$submitter_id, 20))
}

# ------------------------------------------------------------
# 8. Check whether protein sample names resemble
#    patient IDs after common formatting
# ------------------------------------------------------------

protein_ids <- colnames(prad_prot)

cat("\nProtein IDs:\n")
print(protein_ids)

if ("submitter_id" %in% colnames(clinical_prad)) {

  clinical_ids <- clinical_prad$submitter_id

  # Normalize common TCGA punctuation differences
  protein_patient_ids <- gsub(
    "\\.",
    "-",
    protein_ids
  )

  protein_patient_ids <- sub(
    "-[0-9]{2}[A-Z].*$",
    "",
    protein_patient_ids
  )

  clinical_ids_clean <- gsub(
    "\\.",
    "-",
    clinical_ids
  )

  cat("\nNormalized protein IDs:\n")
  print(protein_patient_ids)

  cat("\nNumber of direct patient-ID matches:\n")

  matches <- protein_patient_ids[
    protein_patient_ids %in% clinical_ids_clean
  ]

  print(length(matches))

  cat("\nMatched IDs:\n")
  print(matches)
}


# ------------------------------------------------------------
# 9. Check BRCA protein pairing
# ------------------------------------------------------------

cat("\n================ BRCA PROTEIN PAIRING ================\n")

brca_pairing <- read.csv(
  "01_RawData/Breast/brca_protein_sample_pairing.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Dimensions:\n")
print(dim(brca_pairing))

cat("\nUnique BRCA patients:\n")
print(length(unique(brca_pairing$patient_id)))

cat("\nSample types:\n")
print(table(brca_pairing$sample_type))

cat("\nFirst 20 patient/sample mappings:\n")
print(head(brca_pairing, 20))


# ------------------------------------------------------------
# 10. Save current audit
# ------------------------------------------------------------

write.csv(
  brca_complete,
  "05_FinalResults/Step5_BRCA_complete_omics_patients.csv",
  row.names = FALSE
)

write.csv(
  prad_manifest,
  "05_FinalResults/Step5_PRAD_manifest_audit.csv",
  row.names = FALSE
)

cat("\n============================================================\n")
cat("STEP 5.0G COMPLETE\n")
cat("============================================================\n")