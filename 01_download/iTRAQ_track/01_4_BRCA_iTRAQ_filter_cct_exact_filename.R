# ================================================
# CORRECT PROTEOMICS FILTER
# Using the right file: the .cct file (7.0 MB)
# ================================================

cat("========================================\n")
cat("   LOADING CORRECT PROTEOMICS FILE\n")
cat("========================================\n\n")

# -----------------------------------------------
# STEP 1: Load the correct CCT file
# -----------------------------------------------
cat("STEP 1: Loading CCT file...\n\n")

cct_file <- "01_RawData/Breast/Human__TCGA_BRCA__BI__Proteome__QExact__01_28_2016__BI__Gene__CDAP_iTRAQ_UnsharedLogRatio_r2.cct"

cat(sprintf("File: %s\n", basename(cct_file)))
cat(sprintf("Size: %.1f MB\n\n",
            file.size(cct_file)/1024/1024))

brca_protein_raw <- read.table(
  cct_file,
  header       = TRUE,
  sep          = "\t",
  row.names    = 1,
  check.names  = FALSE,
  comment.char = ""
)

cat(sprintf("✅ Loaded!\n"))
cat(sprintf("  Proteins : %d\n",
            nrow(brca_protein_raw)))
cat(sprintf("  Samples  : %d\n",
            ncol(brca_protein_raw)))

# Show first few sample names
cat("\nSample name examples:\n")
print(head(colnames(brca_protein_raw), 8))

# -----------------------------------------------
# STEP 2: Load your early stage barcodes
# -----------------------------------------------
cat("\nSTEP 2: Loading early stage barcodes...\n\n")

clinical_brca <- read.csv(
  "01_RawData/Breast/clinical_brca.csv")
early_barcodes_brca <- clinical_brca$submitter_id

cat(sprintf("Early-stage patients: %d\n",
            length(early_barcodes_brca)))
cat("Stage distribution:\n")
print(table(clinical_brca$ajcc_pathologic_stage))

# Show barcode examples
cat("\nYour barcode examples:\n")
print(head(early_barcodes_brca, 5))

# -----------------------------------------------
# STEP 3: Standardize and match barcodes
# -----------------------------------------------
cat("\nSTEP 3: Matching barcodes...\n\n")

# Standardize protein column names
# Replace dots with dashes
prot_std  <- gsub("\\.", "-",
                   colnames(brca_protein_raw))

cat("Protein barcode examples (standardized):\n")
print(head(prot_std, 5))

# Match on first 12 characters
# TCGA-XX-XXXX
prot_12 <- substr(prot_std, 1, 12)
clin_12 <- substr(early_barcodes_brca, 1, 12)

cat("\nProtein 12-char examples:\n")
print(head(prot_12, 5))
cat("\nClinical 12-char examples:\n")
print(head(clin_12, 5))

# Count matches
matched_idx <- prot_12 %in% clin_12
n_matched   <- sum(matched_idx)

cat(sprintf("\nTotal protein samples : %d\n",
            ncol(brca_protein_raw)))
cat(sprintf("Early-stage matches   : %d\n",
            n_matched))
cat(sprintf("Percentage matched    : %.1f%%\n",
            100 * n_matched /
            ncol(brca_protein_raw)))

# -----------------------------------------------
# STEP 4: Filter to early stage
# -----------------------------------------------
cat("\nSTEP 4: Filtering...\n\n")

if(n_matched >= 10){

  # Filter columns to early stage only
  brca_protein_early <- brca_protein_raw[,
                          matched_idx]

  cat(sprintf("✅ Filtered successfully!\n"))
  cat(sprintf("  Proteins         : %d\n",
              nrow(brca_protein_early)))
  cat(sprintf("  Early-stage samples: %d\n",
              ncol(brca_protein_early)))

  # Verify stages of matched samples
  matched_barcodes <- prot_std[matched_idx]
  matched_12       <- substr(matched_barcodes,
                              1, 12)
  matched_clinical <- clinical_brca[
    substr(clinical_brca$submitter_id,
           1, 12) %in% matched_12, ]

  cat("\nStage verification of matched samples:\n")
  print(table(
    matched_clinical$ajcc_pathologic_stage))

} else {

  cat(sprintf(
    "⚠️  Only %d barcode matches\n", n_matched))
  cat("Keeping all tumour samples instead\n\n")

  # Tumour = position 14-15 is "01"
  sample_type <- substr(prot_std, 14, 15)
  cat("Sample types found:\n")
  print(table(sample_type))

  tumour_idx         <- sample_type == "01"
  brca_protein_early <- brca_protein_raw[,
                          tumour_idx]

  cat(sprintf("\nTumour samples kept: %d\n",
              ncol(brca_protein_early)))
}

# -----------------------------------------------
# STEP 5: Create cancer vs normal labels
# -----------------------------------------------
cat("\nSTEP 5: Creating group labels...\n\n")

early_std  <- gsub("\\.", "-",
  colnames(brca_protein_early))
early_type <- substr(early_std, 14, 15)

brca_prot_groups <- ifelse(
  early_type == "11",
  "Normal",
  "EarlyBreastCancer"
)

cat("Group distribution:\n")
print(table(brca_prot_groups))

# -----------------------------------------------
# STEP 6: Delete wrong files and save correctly
# -----------------------------------------------
cat("\nSTEP 6: Removing wrong files...\n\n")

# Remove the incorrectly saved files
wrong_files <- c(
  "01_RawData/Breast/brca_protein_early.csv",
  "01_RawData/Breast/brca_protein_raw.csv",
  "01_RawData/Breast/brca_meth_matrix.csv"
)

for(wf in wrong_files){
  if(file.exists(wf)){
    file.remove(wf)
    cat(sprintf("  🗑️  Removed: %s\n",
                basename(wf)))
  }
}

cat("\nSTEP 6: Saving correct files...\n\n")

# Save early stage proteomics
write.csv(
  brca_protein_early,
  "01_RawData/Breast/brca_protein_early.csv"
)
cat(sprintf("✅ SAVED: brca_protein_early.csv  %.1f MB\n",
    file.size("01_RawData/Breast/brca_protein_early.csv")
    /1024/1024))

# Save group labels
saveRDS(
  brca_prot_groups,
  "01_RawData/Breast/brca_prot_groups.rds"
)
cat("✅ SAVED: brca_prot_groups.rds\n")

# -----------------------------------------------
# FINAL SUMMARY
# -----------------------------------------------
cat("\n========================================\n")
cat("   STEP 1.4 COMPLETE SUMMARY\n")
cat("========================================\n\n")

cat(sprintf("Source file    : LinkedOmics CCT\n"))
cat(sprintf("Total proteins : %d\n",
            nrow(brca_protein_early)))
cat(sprintf("Total samples  : %d\n",
            ncol(brca_protein_early)))
cat(sprintf("Cancer samples : %d\n",
  sum(brca_prot_groups == "EarlyBreastCancer")))
cat(sprintf("Normal samples : %d\n",
  sum(brca_prot_groups == "Normal")))
cat(sprintf("File size      : %.1f MB\n",
  file.size("01_RawData/Breast/brca_protein_early.csv")
  /1024/1024))

cat("\n========================================\n")
cat("  ✅ STEP 1.4 COMPLETE\n")
cat("  Continue to Step 1.5\n")
cat("  Breast metabolomics download\n")
cat("========================================\n")