# ================================================
# LOAD AND FILTER PROTEOMICS
# Using index 11 — avoids filename problems
# ================================================

cat("========================================\n")
cat("   PROTEOMICS EARLY STAGE FILTER\n")
cat("========================================\n\n")

# Get all files
all_files <- list.files(
  "01_RawData/Breast/",
  full.names = TRUE
)

# Load file 11 — the CCT proteomics file
cct_file <- all_files[11]

cat(sprintf("Loading: %s\n", basename(cct_file)))
cat(sprintf("Size   : %.1f MB\n\n",
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
cat(sprintf("  Samples  : %d\n\n",
            ncol(brca_protein_raw)))

cat("Sample name examples:\n")
print(head(colnames(brca_protein_raw), 8))

# -----------------------------------------------
# Load early stage barcodes
# -----------------------------------------------
cat("\nLoading early stage barcodes...\n")

clinical_brca <- read.csv(
  "01_RawData/Breast/clinical_brca.csv")
early_barcodes_brca <- clinical_brca$submitter_id

cat(sprintf("Early-stage patients: %d\n\n",
            length(early_barcodes_brca)))

cat("Your barcode examples:\n")
print(head(early_barcodes_brca, 5))

# -----------------------------------------------
# Match barcodes
# -----------------------------------------------
cat("\nMatching barcodes...\n\n")

# Standardize protein names
# Replace dots with dashes
prot_std <- gsub("\\.", "-",
                  colnames(brca_protein_raw))

cat("Protein barcode examples (standardized):\n")
print(head(prot_std, 5))

# Match on 12 characters: TCGA-XX-XXXX
prot_12 <- substr(prot_std, 1, 12)
clin_12 <- substr(early_barcodes_brca, 1, 12)

matched_idx <- prot_12 %in% clin_12
n_matched   <- sum(matched_idx)

cat(sprintf("\nTotal protein samples : %d\n",
            ncol(brca_protein_raw)))
cat(sprintf("Early-stage matches   : %d\n",
            n_matched))
cat(sprintf("Percentage matched    : %.1f%%\n\n",
            100 * n_matched /
            ncol(brca_protein_raw)))

# -----------------------------------------------
# Filter to early stage
# -----------------------------------------------
cat("Filtering to early stage...\n\n")

if(n_matched >= 10){

  brca_protein_early <- brca_protein_raw[,
                          matched_idx]

  cat("✅ Filter successful!\n")
  cat(sprintf("  Proteins           : %d\n",
              nrow(brca_protein_early)))
  cat(sprintf("  Early-stage samples: %d\n\n",
              ncol(brca_protein_early)))

  # Verify stages
  matched_12       <- prot_12[matched_idx]
  matched_clinical <- clinical_brca[
    substr(clinical_brca$submitter_id,
           1, 12) %in% matched_12, ]

  cat("Stage verification:\n")
  print(table(
    matched_clinical$ajcc_pathologic_stage))

} else {

  cat(sprintf("⚠️  Only %d matches\n\n",
              n_matched))
  cat("Keeping all tumour samples (type 01)...\n")

  sample_type <- substr(prot_std, 14, 15)
  cat("Sample types in data:\n")
  print(table(sample_type))

  tumour_idx         <- sample_type == "01"
  brca_protein_early <- brca_protein_raw[,
                          tumour_idx]

  cat(sprintf("\n✅ Tumour samples: %d\n",
              ncol(brca_protein_early)))
}

# -----------------------------------------------
# Create group labels
# -----------------------------------------------
cat("\nCreating group labels...\n\n")

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
# Delete wrong large files
# -----------------------------------------------
cat("\nRemoving incorrectly saved files...\n")

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

# -----------------------------------------------
# Save correctly
# -----------------------------------------------
cat("\nSaving correct files...\n\n")

write.csv(
  brca_protein_early,
  "01_RawData/Breast/brca_protein_early.csv"
)

size_saved <- round(
  file.size(
    "01_RawData/Breast/brca_protein_early.csv"
  ) / 1024 / 1024, 1
)

cat(sprintf("✅ SAVED: brca_protein_early.csv\n"))
cat(sprintf("   Size : %.1f MB\n", size_saved))

saveRDS(
  brca_prot_groups,
  "01_RawData/Breast/brca_prot_groups.rds"
)
cat("✅ SAVED: brca_prot_groups.rds\n")

# -----------------------------------------------
# Final summary
# -----------------------------------------------
cat("\n========================================\n")
cat("   STEP 1.4 SUMMARY\n")
cat("========================================\n\n")

cat(sprintf("Source  : LinkedOmics CCT file\n"))
cat(sprintf("Proteins: %d\n",
            nrow(brca_protein_early)))
cat(sprintf("Samples : %d\n",
            ncol(brca_protein_early)))
cat(sprintf("Cancer  : %d\n",
  sum(brca_prot_groups == "EarlyBreastCancer")))
cat(sprintf("Normal  : %d\n",
  sum(brca_prot_groups == "Normal")))
cat(sprintf("File    : %.1f MB\n", size_saved))

# Size check
if(size_saved < 100){
  cat("\n✅ File size looks correct\n")
  cat("   (small file = right dataset)\n")
} else {
  cat("\n⚠️  File is large — check content\n")
}

cat("\n========================================\n")
cat("  ✅ STEP 1.4 COMPLETE\n")
cat("  Continue to Step 1.5\n")
cat("  Breast metabolomics\n")
cat("========================================\n")