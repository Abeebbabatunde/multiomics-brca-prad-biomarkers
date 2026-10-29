# ============================================================
# STEP 2.3: Breast Proteomics "Normalization" 
# (QC + missing value handling on pre-normalized RPPA log-ratios)
# ============================================================

library(dplyr)
library(tidyr)

brca_protein_raw <- read.csv("01_RawData/Breast/brca_protein_raw.csv", stringsAsFactors = FALSE)

# ---- 1. Quantify missingness per protein and per sample ----
cat("STEP 1: Missingness audit\n")

missing_by_protein <- brca_protein_raw %>%
  group_by(peptide_target) %>%
  summarise(n_missing = sum(is.na(value)), pct_missing = mean(is.na(value)) * 100) %>%
  arrange(desc(pct_missing))

missing_by_sample <- brca_protein_raw %>%
  group_by(sample_barcode, sample_type) %>%
  summarise(n_missing = sum(is.na(value)), pct_missing = mean(is.na(value)) * 100) %>%
  arrange(desc(pct_missing))

cat("Proteins with most missingness:\n")
print(head(missing_by_protein, 10))

cat("\nSamples with most missingness:\n")
print(head(missing_by_sample, 10))

# ---- 2. Decide a missingness threshold and filter proteins ----
# Common RPPA practice: drop proteins missing in >20% of samples
threshold_pct <- 20

proteins_to_keep <- missing_by_protein$peptide_target[missing_by_protein$pct_missing <= threshold_pct]
proteins_dropped <- missing_by_protein$peptide_target[missing_by_protein$pct_missing > threshold_pct]

cat("\nProteins kept (≤", threshold_pct, "% missing):", length(proteins_to_keep), "\n")
cat("Proteins dropped (>", threshold_pct, "% missing):", length(proteins_dropped), "\n")
if (length(proteins_dropped) > 0) print(proteins_dropped)

brca_protein_filtered <- brca_protein_raw %>%
  filter(peptide_target %in% proteins_to_keep)

# ---- 3. Impute remaining sparse NAs (small amount) using per-protein median ----
cat("\nSTEP 3: Imputing remaining missing values (per-protein median)\n")

brca_protein_normalized <- brca_protein_filtered %>%
  group_by(peptide_target) %>%
  mutate(value = ifelse(is.na(value), median(value, na.rm = TRUE), value)) %>%
  ungroup()

cat("Remaining NAs after imputation:", sum(is.na(brca_protein_normalized$value)), "\n")

# ---- 4. Confirm tumor/normal pairing is still intact post-filtering ----
cat("\nSTEP 4: Re-verify pairing intact\n")
table(brca_protein_normalized$sample_type)
n_distinct_patients <- length(unique(substr(brca_protein_normalized$sample_barcode, 1, 12)))
cat("Distinct patients remaining:", n_distinct_patients, "(should still be 30)\n")

# ---- 5. Save Step 2.3 output ----
norm_dir <- "02_Normalized/Breast"  # adjust if your pipeline uses a different folder name
dir.create(norm_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(brca_protein_normalized, 
          file.path(norm_dir, "brca_protein_normalized.csv"), 
          row.names = FALSE)

cat("\n✅ Saved:", file.path(norm_dir, "brca_protein_normalized.csv"), "\n")
cat("Final dimensions:", nrow(brca_protein_normalized), "rows,", 
    length(unique(brca_protein_normalized$peptide_target)), "proteins,",
    length(unique(brca_protein_normalized$sample_barcode)), "samples\n")
