# ============================================================
# STEP 5.0H — AUDIT ACTUAL BRCA THREE-OMICS SAMPLE IDs
# ============================================================

cat("\n============================================================\n")
cat("STEP 5.0H — BRCA THREE-OMICS MATRIX AUDIT\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Load the 31 harmonized BRCA patients
# ------------------------------------------------------------

brca_complete <- read.csv(
  "05_FinalResults/Step5_BRCA_complete_omics_patients.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\nNumber of harmonized BRCA patients:",
    nrow(brca_complete), "\n")

cat("\nHarmonized patient IDs:\n")
print(brca_complete$patient_id)


# ------------------------------------------------------------
# 2. BRCA RNA
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("BRCA RNA\n")
cat("============================================================\n")

brca_rna <- readRDS(
  "01_RawData/Breast/brca_rna.rds"
)

cat("Class:", class(brca_rna), "\n")
cat("Dimensions:", paste(dim(brca_rna), collapse = " x "), "\n")

rna_ids <- colnames(brca_rna)

cat("\nFirst 10 RNA sample IDs:\n")
print(head(rna_ids, 10))


# ------------------------------------------------------------
# 3. Extract patient IDs from RNA sample IDs
# ------------------------------------------------------------

rna_patient_ids <- sub(
  "-[0-9]{2}[A-Z].*$",
  "",
  rna_ids
)

rna_patient_ids <- sub(
  "\\.[0-9]+$",
  "",
  rna_patient_ids
)

cat("\nNumber of unique RNA patient IDs:",
    length(unique(rna_patient_ids)), "\n")


# ------------------------------------------------------------
# 4. BRCA methylation raw object
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("BRCA METHYLATION\n")
cat("============================================================\n")

brca_meth <- readRDS(
  "01_RawData/Breast/brca_meth.rds"
)

cat("Class:", class(brca_meth), "\n")
cat("Dimensions:", paste(dim(brca_meth), collapse = " x "), "\n")

meth_ids <- colnames(brca_meth)

cat("\nFirst 10 methylation sample IDs:\n")
print(head(meth_ids, 10))


# ------------------------------------------------------------
# 5. Extract patient IDs from methylation
# ------------------------------------------------------------

meth_patient_ids <- sub(
  "-[0-9]{2}[A-Z].*$",
  "",
  meth_ids
)

meth_patient_ids <- sub(
  "\\.[0-9]+$",
  "",
  meth_patient_ids
)

cat("\nNumber of unique methylation patient IDs:",
    length(unique(meth_patient_ids)), "\n")


# ------------------------------------------------------------
# 6. Locate BRCA proteomics files
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("BRCA PROTEOMICS FILE SEARCH\n")
cat("============================================================\n")

brca_files <- list.files(
  "01_RawData/Breast",
  recursive = TRUE,
  full.names = TRUE
)

brca_prot_files <- brca_files[
  grepl(
    "prot|protein|proteom",
    basename(brca_files),
    ignore.case = TRUE
  )
]

cat("\nPotential BRCA proteomics files:\n")
print(brca_prot_files)


# ------------------------------------------------------------
# 7. Inspect BRCA protein pairing
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("BRCA PROTEIN SAMPLE PAIRING\n")
cat("============================================================\n")

brca_pairing <- read.csv(
  "01_RawData/Breast/brca_protein_sample_pairing.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Dimensions:",
    paste(dim(brca_pairing), collapse = " x "),
    "\n")

print(head(brca_pairing, 20))

cat("\nUnique patients:",
    length(unique(brca_pairing$patient_id)),
    "\n")

cat("\nUnique sample barcodes:",
    length(unique(brca_pairing$sample_barcode)),
    "\n")

cat("\nSample types:\n")
print(table(brca_pairing$sample_type))


# ------------------------------------------------------------
# 8. Compare the 31 patients across RNA and methylation
# ------------------------------------------------------------

target_ids <- unique(brca_complete$patient_id)

rna_match <- target_ids %in% unique(rna_patient_ids)

meth_match <- target_ids %in% unique(meth_patient_ids)

prot_match <- target_ids %in%
  unique(brca_pairing$patient_id)


cat("\n============================================================\n")
cat("31-PATIENT MATCH RESULTS\n")
cat("============================================================\n")

cat("\nRNA matches:",
    sum(rna_match),
    "of",
    length(target_ids),
    "\n")

cat("Methylation matches:",
    sum(meth_match),
    "of",
    length(target_ids),
    "\n")

cat("Proteomics matches:",
    sum(prot_match),
    "of",
    length(target_ids),
    "\n")


# ------------------------------------------------------------
# 9. Detailed matching table
# ------------------------------------------------------------

brca_id_audit <- data.frame(
  patient_id = target_ids,
  RNA = rna_match,
  Methylation = meth_match,
  Proteomics = prot_match,
  stringsAsFactors = FALSE
)

cat("\nDetailed audit:\n")
print(brca_id_audit)


# ------------------------------------------------------------
# 10. Identify any missing IDs
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("MISSING IDs\n")
cat("============================================================\n")

cat("\nMissing from RNA:\n")
print(
  brca_id_audit$patient_id[
    !brca_id_audit$RNA
  ]
)

cat("\nMissing from methylation:\n")
print(
  brca_id_audit$patient_id[
    !brca_id_audit$Methylation
  ]
)

cat("\nMissing from proteomics:\n")
print(
  brca_id_audit$patient_id[
    !brca_id_audit$Proteomics
  ]
)


# ------------------------------------------------------------
# 11. Save audit
# ------------------------------------------------------------

write.csv(
  brca_id_audit,
  "05_FinalResults/Step5_BRCA_three_omics_ID_audit.csv",
  row.names = FALSE
)


cat("\n============================================================\n")
cat("STEP 5.0H COMPLETE\n")
cat("Audit saved to:\n")
cat("05_FinalResults/Step5_BRCA_three_omics_ID_audit.csv\n")
cat("============================================================\n")