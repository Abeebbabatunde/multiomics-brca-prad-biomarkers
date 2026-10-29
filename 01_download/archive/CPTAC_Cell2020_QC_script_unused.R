# ============================================================
# Step 2.3 - QC Check: CPTAC Breast Cancer Proteomics Data
# Folder: MyThesis_MultiOmics/01_RawData/Breast/
#         breast cancer (CPTAC, Cell 2020)
# Purpose: Verify file integrity + confirm presence of
#          early-stage tumor AND normal samples
# ============================================================

base_dir <- "MyThesis_MultiOmics/01_RawData/Breast/breast cancer (CPTAC, Cell 2020)"
# ^ adjust capitalization/spacing to match exactly what you see in your file explorer

cat("STEP 1: Folder contents\n")
all_files <- list.files(base_dir, recursive = TRUE, full.names = TRUE)
print(all_files)

if (length(all_files) == 0) stop("No files found - check base_dir path is correct")

# STEP 2: Locate key files by pattern (cBioPortal-style naming)
find_file <- function(pattern) {
  hit <- all_files[grepl(pattern, basename(all_files), ignore.case = TRUE)]
  if (length(hit) == 0) return(NA)
  hit[1]
}

f_clinical_sample  <- find_file("clinical_sample")
f_clinical_patient <- find_file("clinical_patient")
f_protein          <- find_file("protein_quantification|protein_expression|data_protein")

cat("\nclinical_sample file: ", f_clinical_sample, "\n")
cat("clinical_patient file:", f_clinical_patient, "\n")
cat("protein data file:    ", f_protein, "\n")

# STEP 3: Basic file integrity checks (not empty, not corrupted)
check_file <- function(path, label) {
  if (is.na(path)) { cat(label, "- MISSING\n"); return(invisible(NULL)) }
  sz <- file.size(path)
  cat(label, "- found,", round(sz/1024, 1), "KB\n")
  if (sz == 0) warning(paste(label, "is empty (0 bytes) - re-download"))
}
check_file(f_clinical_sample, "Clinical sample file")
check_file(f_clinical_patient, "Clinical patient file")
check_file(f_protein, "Protein quantification file")

# STEP 4: Read clinical sample data, inspect sample type & stage
clin_sample <- NULL
if (!is.na(f_clinical_sample)) {
  clin_sample <- read.delim(f_clinical_sample, comment.char = "#",
                             stringsAsFactors = FALSE, check.names = FALSE)
  cat("\nClinical sample columns:\n")
  print(colnames(clin_sample))

  type_col <- grep("SAMPLE_TYPE|TISSUE|SPECIMEN", colnames(clin_sample),
                    ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(type_col)) {
    cat("\nSample type counts (column:", type_col, ") - look for Tumor AND Normal:\n")
    print(table(clin_sample[[type_col]]))
  } else {
    cat("\nNo obvious sample-type column auto-detected. Inspect manually:\n")
    print(head(clin_sample, 3))
  }

  stage_col <- grep("STAGE", colnames(clin_sample), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(stage_col)) {
    cat("\nTumor stage counts (column:", stage_col, ") - look for Stage I/II:\n")
    print(table(clin_sample[[stage_col]]))
  }
}

# STEP 5: Check patient-level clinical file for stage if not found above
if (!is.na(f_clinical_patient)) {
  clin_patient <- read.delim(f_clinical_patient, comment.char = "#",
                              stringsAsFactors = FALSE, check.names = FALSE)
  stage_col_p <- grep("STAGE", colnames(clin_patient), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(stage_col_p)) {
    cat("\nPatient-level stage counts (column:", stage_col_p, "):\n")
    print(table(clin_patient[[stage_col_p]]))
  }
}

# STEP 6: Read protein matrix - dimensions, missingness, sample ID overlap
if (!is.na(f_protein)) {
  protein_data <- read.delim(f_protein, stringsAsFactors = FALSE, check.names = FALSE)
  cat("\nProtein data dimensions:", nrow(protein_data), "rows x", ncol(protein_data), "cols\n")
  cat("First 10 column names:\n")
  print(head(colnames(protein_data), 10))

  numeric_cols <- sapply(protein_data, is.numeric)
  pct_na <- mean(is.na(protein_data[, numeric_cols]))
  cat("Overall missing value rate in protein matrix:", round(pct_na * 100, 2), "%\n")

  if (!is.null(clin_sample) && "SAMPLE_ID" %in% colnames(clin_sample)) {
    overlap <- intersect(colnames(protein_data), clin_sample$SAMPLE_ID)
    cat("Protein columns matching clinical SAMPLE_ID:", length(overlap),
        "out of", ncol(protein_data), "columns\n")
  }
}

cat("\n=== SUMMARY - what to check ===\n")
cat("1. Sample-type table above must show BOTH a tumor category and a normal category.\n")
cat("2. Stage table must include Stage I and/or Stage II entries (your early-detection criterion).\n")
cat("3. Missing-value rate above ~30-40% may need imputation or filtering before proceeding.\n")
cat("4. If sample-type/stage columns weren't auto-detected, check the printed column names\n")
cat("   and column previews above, then re-run table() on the correct column manually.\n")