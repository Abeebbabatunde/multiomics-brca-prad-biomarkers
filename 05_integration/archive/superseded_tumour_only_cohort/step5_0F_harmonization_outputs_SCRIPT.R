# ============================================================
# STEP 5.0F — INSPECT CLINICAL HARMONIZATION OUTPUTS
# ============================================================

cat("\n============================================================\n")
cat("1. HARMONIZED SAMPLE MANIFEST\n")
cat("============================================================\n")

manifest <- read.csv(
  "00.5_clinical_metadata_harmonization/sample_manifest_harmonized.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Dimensions:\n")
print(dim(manifest))

cat("\nColumn names:\n")
print(colnames(manifest))

cat("\nFirst 10 rows:\n")
print(head(manifest, 10))


# ============================================================
# 2. ATTRITION REPORT\n# ============================================================

cat("\n============================================================\n")
cat("2. ATTRITION REPORT\n")
cat("============================================================\n")

attrition <- read.csv(
  "00.5_clinical_metadata_harmonization/attrition_report.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Dimensions:\n")
print(dim(attrition))

cat("\nColumn names:\n")
print(colnames(attrition))

cat("\nFull attrition report:\n")
print(attrition)


# ============================================================
# 3. BRCA PROTEIN SAMPLE PAIRING\n# ============================================================

cat("\n============================================================\n")
cat("3. BRCA PROTEIN SAMPLE PAIRING\n")
cat("============================================================\n")

brca_pairing <- read.csv(
  "01_RawData/Breast/brca_protein_sample_pairing.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Dimensions:\n")
print(dim(brca_pairing))

cat("\nColumn names:\n")
print(colnames(brca_pairing))

cat("\nFirst 20 rows:\n")
print(head(brca_pairing, 20))


# ============================================================
# 4. BRCA PROTEOMICS METADATA\n# ============================================================

cat("\n============================================================\n")
cat("4. BRCA PROTEOMICS METADATA\n")
cat("============================================================\n")

brca_prot_meta_csv <- read.csv(
  "01_RawData/Breast/brca_prot_meta_full.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Dimensions:\n")
print(dim(brca_prot_meta_csv))

cat("\nColumn names:\n")
print(colnames(brca_prot_meta_csv))

cat("\nFirst 10 rows:\n")
print(head(brca_prot_meta_csv, 10))


# ============================================================
# 5. BRCA CLINICAL DATA\n# ============================================================

cat("\n============================================================\n")
cat("5. BRCA CLINICAL DATA\n")
cat("============================================================\n")

clinical_brca <- read.csv(
  "01_RawData/Breast/clinical_brca.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Dimensions:\n")
print(dim(clinical_brca))

cat("\nColumn names:\n")
print(colnames(clinical_brca))

cat("\nFirst 5 rows:\n")
print(head(clinical_brca, 5))


# ============================================================
# 6. PRAD CLINICAL DATA\n# ============================================================

cat("\n============================================================\n")
cat("6. PRAD CLINICAL DATA\n")
cat("============================================================\n")

clinical_prad <- read.csv(
  "01_RawData/Prostate/clinical_prad.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Dimensions:\n")
print(dim(clinical_prad))

cat("\nColumn names:\n")
print(colnames(clinical_prad))

cat("\nFirst 5 rows:\n")
print(head(clinical_prad, 5))


# ============================================================
# 7. FIND THE 00.5 HARMONIZATION SCRIPT(S)
# ============================================================

cat("\n============================================================\n")
cat("7. 00.5 HARMONIZATION SCRIPTS\n")
cat("============================================================\n")

harmonization_files <- list.files(
  "00.5_clinical_metadata_harmonization",
  recursive = TRUE,
  full.names = TRUE
)

print(harmonization_files)