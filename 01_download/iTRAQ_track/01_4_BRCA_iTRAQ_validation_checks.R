# ================================================
# VALIDATION OF EARLY-STAGE PROTEOMICS
# ================================================

cat("\n========================================\n")
cat("   VALIDATING DATASET\n")
cat("========================================\n\n")

# -----------------------------------------------
# CHECK 1 — Matrix dimensions
# -----------------------------------------------
cat("CHECK 1: Matrix dimensions\n")

dims <- dim(brca_protein_early)

cat(sprintf("Rows    (proteins): %d\n", dims[1]))
cat(sprintf("Columns (samples) : %d\n\n", dims[2]))

if(dims[1] < 100){
  cat("❌ Too few proteins detected\n")
} else {
  cat("✅ Protein count looks realistic\n")
}

if(dims[2] < 10){
  cat("❌ Too few samples detected\n")
} else {
  cat("✅ Sample count looks realistic\n")
}

# -----------------------------------------------
# CHECK 2 — File size
# -----------------------------------------------
cat("\nCHECK 2: File size\n")

f <- "01_RawData/Breast/brca_protein_early.csv"

size_mb <- file.size(f)/1024/1024

cat(sprintf("File size: %.2f MB\n", size_mb))

if(size_mb < 1){
  cat("❌ File suspiciously small\n")
} else {
  cat("✅ File size looks normal\n")
}

# -----------------------------------------------
# CHECK 3 — Preview data
# -----------------------------------------------
cat("\nCHECK 3: Data preview\n")

print(
  brca_protein_early[
    1:min(5,nrow(brca_protein_early)),
    1:min(5,ncol(brca_protein_early))
  ]
)

# -----------------------------------------------
# CHECK 4 — Missing values
# -----------------------------------------------
cat("\nCHECK 4: Missing values\n")

na_count <- sum(is.na(brca_protein_early))

cat(sprintf("Missing values: %d\n", na_count))

na_percent <- 100 * na_count /
  (nrow(brca_protein_early) *
   ncol(brca_protein_early))

cat(sprintf("Missing percentage: %.2f%%\n",
            na_percent))

if(na_percent > 50){
  cat("⚠️ Very high missingness\n")
} else {
  cat("✅ Missingness acceptable\n")
}

# -----------------------------------------------
# CHECK 5 — Numeric values
# -----------------------------------------------
cat("\nCHECK 5: Numeric validation\n")

numeric_check <- all(
  sapply(brca_protein_early, is.numeric)
)

cat(sprintf("All numeric: %s\n",
            numeric_check))

if(numeric_check){
  cat("✅ Protein values are numeric\n")
} else {
  cat("⚠️ Some columns are not numeric\n")
}

# -----------------------------------------------
# CHECK 6 — Expression distribution
# -----------------------------------------------
cat("\nCHECK 6: Expression ranges\n")

vals <- as.numeric(as.matrix(brca_protein_early))

vals <- vals[!is.na(vals)]

cat(sprintf("Minimum value : %.3f\n", min(vals)))
cat(sprintf("Maximum value : %.3f\n", max(vals)))
cat(sprintf("Mean value    : %.3f\n", mean(vals)))

if(max(vals) == 0){
  cat("❌ All values are zero\n")
} else {
  cat("✅ Expression values detected\n")
}

# -----------------------------------------------
# CHECK 7 — Barcode validation
# -----------------------------------------------
cat("\nCHECK 7: Barcode validation\n")

sample_examples <- head(
  colnames(brca_protein_early), 10)

print(sample_examples)

valid_tcga <- grepl("^TCGA",
                    sample_examples)

cat(sprintf(
  "Valid TCGA barcodes: %d/%d\n",
  sum(valid_tcga),
  length(valid_tcga)
))

# -----------------------------------------------
# CHECK 8 — Stage verification
# -----------------------------------------------
cat("\nCHECK 8: Stage verification\n")

if(exists("matched_clinical")){

  print(
    table(
      matched_clinical$ajcc_pathologic_stage
    )
  )

  cat("\n✅ Stage matching exists\n")

} else {

  cat("⚠️ matched_clinical object missing\n")
  cat("Stage verification unavailable\n")
}

# -----------------------------------------------
# CHECK 9 — Empty rows/columns
# -----------------------------------------------
cat("\nCHECK 9: Empty rows/columns\n")

empty_rows <- sum(
  rowSums(is.na(brca_protein_early)) ==
    ncol(brca_protein_early)
)

empty_cols <- sum(
  colSums(is.na(brca_protein_early)) ==
    nrow(brca_protein_early)
)

cat(sprintf("Empty proteins : %d\n",
            empty_rows))
cat(sprintf("Empty samples  : %d\n",
            empty_cols))

if(empty_rows == 0 & empty_cols == 0){
  cat("✅ No completely empty rows/columns\n")
}

# -----------------------------------------------
# FINAL DECISION
# -----------------------------------------------
cat("\n========================================\n")
cat("   FINAL VALIDATION\n")
cat("========================================\n\n")

if(
  dims[1] > 100 &
  dims[2] >= 10 &
  size_mb > 1 &
  max(vals) > 0
){

  cat("✅ DATASET APPEARS VALID\n")
  cat("✅ Proteomics successfully loaded\n")
  cat("✅ Early-stage filtering likely worked\n")
  cat("✅ Safe to continue analysis\n")

} else {

  cat("⚠️ DATASET MAY HAVE ISSUES\n")
  cat("Inspect the checks above carefully\n")
}

cat("\n========================================\n")