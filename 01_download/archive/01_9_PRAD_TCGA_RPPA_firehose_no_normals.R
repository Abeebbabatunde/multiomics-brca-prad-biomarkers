# ================================================
# LOAD PROSTATE RPPA PROTEOMICS
# ================================================

cat("========================================\n")
cat("   STEP 1.9: PROSTATE PROTEOMICS\n")
cat("========================================\n\n")

# -----------------------------------------------
# STEP 1: Find the CCT file
# -----------------------------------------------
cat("STEP 1: Finding CCT file...\n\n")

all_prad_files <- list.files(
  "01_RawData/Prostate/",
  full.names = TRUE
)

# Find CCT file specifically
cct_idx <- grep(
  "\\.cct$|RPPA|Proteome|proteom",
  basename(all_prad_files),
  ignore.case = TRUE
)

cat("Files found:\n")
for(i in seq_along(all_prad_files)){
  cat(sprintf("  [%d] %s  (%.2f MB)\n",
              i,
              basename(all_prad_files[i]),
              file.size(all_prad_files[i])/
              1024/1024))
}

cat(sprintf("\nCCT file index: %s\n\n",
            paste(cct_idx, collapse=", ")))

# Load using index
prot_file <- all_prad_files[cct_idx[1]]

cat(sprintf("Loading: %s\n",
            basename(prot_file)))
cat(sprintf("Size   : %.2f MB\n\n",
            file.size(prot_file)/1024/1024))

# -----------------------------------------------
# STEP 2: Peek at structure first
# -----------------------------------------------
cat("STEP 2: Peeking at file structure...\n\n")

raw_lines <- readLines(prot_file, n = 3)

cat("Line 1 (first 300 chars):\n")
cat(substr(raw_lines[1], 1, 300), "\n\n")

cat("Line 2 (first 300 chars):\n")
cat(substr(raw_lines[2], 1, 300), "\n\n")

n_cols <- length(
  strsplit(raw_lines[1], "\t")[[1]])
cat(sprintf("Columns detected: %d\n\n", n_cols))

# -----------------------------------------------
# STEP 3: Load the file
# -----------------------------------------------
cat("STEP 3: Loading file...\n\n")

prad_protein_raw <- read.table(
  prot_file,
  header       = TRUE,
  sep          = "\t",
  row.names    = 1,
  check.names  = FALSE,
  comment.char = "",
  fill         = TRUE,
  quote        = ""
)

cat(sprintf("✅ Loaded!\n"))
cat(sprintf("  Proteins : %d\n",
            nrow(prad_protein_raw)))
cat(sprintf("  Samples  : %d\n\n",
            ncol(prad_protein_raw)))

cat("Sample name examples:\n")
print(head(colnames(prad_protein_raw), 8))

cat("\nProtein name examples:\n")
print(head(rownames(prad_protein_raw), 8))

cat("\nValue range:\n")
print(range(prad_protein_raw, na.rm=TRUE))

# -----------------------------------------------
# STEP 4: Filter to early stage
# -----------------------------------------------
cat("\nSTEP 4: Filtering to early stage...\n\n")

# Load prostate clinical data
if(!exists("clinical_prad")){
  clinical_prad <- read.csv(
    "01_RawData/Prostate/clinical_prad.csv")
}

early_barcodes_prad <- clinical_prad$submitter_id

cat(sprintf("Early-stage patients: %d\n\n",
            length(early_barcodes_prad)))

# Standardize column names
prot_std <- gsub("\\.", "-",
                  colnames(prad_protein_raw))

# Match on 12 characters
prot_12 <- substr(prot_std, 1, 12)
clin_12 <- substr(early_barcodes_prad, 1, 12)

matched_idx <- prot_12 %in% clin_12
n_matched   <- sum(matched_idx)

cat(sprintf("Total protein samples : %d\n",
            ncol(prad_protein_raw)))
cat(sprintf("Early-stage matches   : %d\n",
            n_matched))
cat(sprintf("Percentage matched    : %.1f%%\n\n",
            100 * n_matched /
            ncol(prad_protein_raw)))

if(n_matched >= 10){

  prad_protein_early <- prad_protein_raw[,
                          matched_idx]
  cat("✅ Filtered to early stage!\n")
  cat(sprintf("  Proteins           : %d\n",
              nrow(prad_protein_early)))
  cat(sprintf("  Early-stage samples: %d\n\n",
              ncol(prad_protein_early)))

  # Verify stages
  matched_12       <- prot_12[matched_idx]
  matched_clinical <- clinical_prad[
    substr(clinical_prad$submitter_id,
           1, 12) %in% matched_12, ]

  cat("Stage verification:\n")
  print(table(
    matched_clinical$ajcc_pathologic_t))

} else {

  cat(sprintf(
    "⚠️  Only %d matches found\n\n", n_matched))
  cat("Keeping all tumour samples...\n")

  # Keep tumour samples using barcode type
  sample_type <- substr(prot_std, 14, 15)
  cat("Sample types:\n")
  print(table(sample_type))

  tumour_idx         <- sample_type == "01"
  prad_protein_early <- prad_protein_raw[,
                          tumour_idx]

  cat(sprintf("\nTumour samples: %d\n",
              ncol(prad_protein_early)))
}

# -----------------------------------------------
# STEP 5: Create group labels
# -----------------------------------------------
cat("\nSTEP 5: Creating group labels...\n\n")

early_std  <- gsub("\\.", "-",
  colnames(prad_protein_early))
early_type <- substr(early_std, 14, 15)

prad_prot_groups <- ifelse(
  early_type == "11",
  "Normal",
  "EarlyProstateCancer"
)

cat("Group distribution:\n")
print(table(prad_prot_groups))

# -----------------------------------------------
# STEP 6: Save
# -----------------------------------------------
cat("\nSTEP 6: Saving...\n\n")

write.csv(
  prad_protein_early,
  "01_RawData/Prostate/prad_protein_early.csv"
)

saveRDS(
  prad_prot_groups,
  "01_RawData/Prostate/prad_prot_groups.rds"
)

size1 <- round(file.size(
  "01_RawData/Prostate/prad_protein_early.csv"
) / 1024 / 1024, 2)

cat(sprintf(
  "✅ SAVED: prad_protein_early.csv  %.2f MB\n",
  size1))
cat("✅ SAVED: prad_prot_groups.rds\n")

# -----------------------------------------------
# FINAL SUMMARY
# -----------------------------------------------
cat("\n========================================\n")
cat("   STEP 1.9 SUMMARY\n")
cat("========================================\n\n")

cat(sprintf("Source   : TCGA-PRAD RPPA\n"))
cat(sprintf("Proteins : %d\n",
            nrow(prad_protein_early)))
cat(sprintf("Samples  : %d\n",
            ncol(prad_protein_early)))
cat(sprintf("Cancer   : %d\n",
  sum(prad_prot_groups ==
      "EarlyProstateCancer")))
cat(sprintf("Normal   : %d\n",
  sum(prad_prot_groups == "Normal")))
cat(sprintf("File size: %.2f MB\n", size1))

cat("\n========================================\n")
cat("  ✅ STEP 1.9 COMPLETE\n")
cat("  Continue to Step 1.10\n")
cat("  Prostate metabolomics\n")
cat("========================================\n")