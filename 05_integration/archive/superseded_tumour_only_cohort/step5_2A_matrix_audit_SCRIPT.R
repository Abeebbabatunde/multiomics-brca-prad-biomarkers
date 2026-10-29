# ============================================================
# STEP 5.2A — AUDIT ACTUAL OMICS MATRICES
# ============================================================

cat("\n============================================================\n")
cat("STEP 5.2A — ACTUAL OMICS MATRIX AUDIT\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Load final cohorts
# ------------------------------------------------------------

brca_3omics <- readRDS(
  "05_FinalResults/Step5_BRCA_3omics_cohort.rds"
)

prad_2omics <- readRDS(
  "05_FinalResults/Step5_PRAD_2omics_cohort.rds"
)

brca_ids <- brca_3omics$patient_id
prad_ids <- prad_2omics$patient_id


# ============================================================
# PART A — BRCA RNA
# ============================================================

cat("\n============================================================\n")
cat("BRCA RNA\n")
cat("============================================================\n")

brca_rna_obj <- readRDS(
  "01_RawData/Breast/brca_rna.rds"
)

cat("Object class:\n")
print(class(brca_rna_obj))

cat("\nDimensions:\n")
print(dim(brca_rna_obj))

cat("\nColumn/sample names — first 10:\n")
print(head(colnames(brca_rna_obj), 10))


# ============================================================
# PART B — BRCA METHYLATION RAW
# ============================================================

cat("\n============================================================\n")
cat("BRCA METHYLATION RAW\n")
cat("============================================================\n")

brca_meth_obj <- readRDS(
  "01_RawData/Breast/brca_meth.rds"
)

cat("Object class:\n")
print(class(brca_meth_obj))

cat("\nDimensions:\n")
print(dim(brca_meth_obj))

cat("\nColumn/sample names — first 10:\n")
print(head(colnames(brca_meth_obj), 10))


# ============================================================
# PART C — BRCA METHYLATION PROCESSED
# ============================================================

cat("\n============================================================\n")
cat("BRCA METHYLATION PROCESSED\n")
cat("============================================================\n")

brca_meth_filtered <- readRDS(
  "02_ProcessedData/Breast/brca_meth_filtered_matrix.rds"
)

cat("Object class:\n")
print(class(brca_meth_filtered))

cat("\nDimensions:\n")
print(dim(brca_meth_filtered))

cat("\nColumn/sample names — first 10:\n")
print(head(colnames(brca_meth_filtered), 10))


# ============================================================
# PART D — BRCA PROTEOMICS FILE SEARCH
# ============================================================

cat("\n============================================================\n")
cat("BRCA PROTEOMICS FILES\n")
cat("============================================================\n")

brca_files <- list.files(
  "01_RawData/Breast",
  recursive = TRUE,
  full.names = TRUE
)

print(brca_files)


# ============================================================
# PART E — PRAD RNA
# ============================================================

cat("\n============================================================\n")
cat("PRAD RNA\n")
cat("============================================================\n")

prad_rna_obj <- readRDS(
  "01_RawData/Prostate/prad_rna.rds"
)

cat("Object class:\n")
print(class(prad_rna_obj))

cat("\nDimensions:\n")
print(dim(prad_rna_obj))

cat("\nColumn/sample names — first 10:\n")
print(head(colnames(prad_rna_obj), 10))


# ============================================================
# PART F — PRAD METHYLATION
# ============================================================

cat("\n============================================================\n")
cat("PRAD METHYLATION\n")
cat("============================================================\n")

prad_meth_obj <- readRDS(
  "01_RawData/Prostate/prad_meth.rds"
)

cat("Object class:\n")
print(class(prad_meth_obj))

cat("\nDimensions:\n")
print(dim(prad_meth_obj))

cat("\nColumn/sample names — first 10:\n")
print(head(colnames(prad_meth_obj), 10))


# ============================================================
# FINAL — COMPARE TARGET IDs WITH RAW OBJECT IDs
# ============================================================

cat("\n============================================================\n")
cat("TARGET PATIENT ID PREVIEW\n")
cat("============================================================\n")

cat("\nBRCA target IDs:\n")
print(head(brca_ids, 10))

cat("\nPRAD target IDs:\n")
print(head(prad_ids, 10))

cat("\n============================================================\n")
cat("STEP 5.2A COMPLETE\n")
cat("============================================================\n")