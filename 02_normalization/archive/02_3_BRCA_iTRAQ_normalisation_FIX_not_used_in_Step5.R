# ================================================
# FIX STEP 2.3 — BREAST PROTEOMICS NORMALIZATION
# ================================================

cat("========================================\n")
cat("   FIXING BREAST PROTEOMICS NORM\n")
cat("========================================\n\n")

# -----------------------------------------------
# STEP 1: Reload protein data correctly
# -----------------------------------------------
cat("STEP 1: Reloading protein data...\n\n")

brca_protein_raw <- read.csv(
  "01_RawData/Breast/brca_protein_early.csv",
  row.names   = 1,
  check.names = FALSE
)

cat(sprintf("Proteins: %d\n", nrow(brca_protein_raw)))
cat(sprintf("Samples : %d\n\n", ncol(brca_protein_raw)))

cat("Sample name examples:\n")
print(head(colnames(brca_protein_raw), 6))

# -----------------------------------------------
# STEP 2: Create group labels from sample names
# LinkedOmics uses TCGA barcodes
# Position 14-15 = sample type
# 01 = tumour, 11 = normal
# -----------------------------------------------
cat("\nSTEP 2: Creating group labels...\n\n")

sample_names <- colnames(brca_protein_raw)
std_names    <- gsub("\\.", "-", sample_names)
sample_type  <- substr(std_names, 14, 15)

cat("Sample types found:\n")
print(table(sample_type))

brca_prot_groups <- ifelse(
  sample_type == "11",
  "Normal",
  "EarlyBreastCancer"
)

cat("\nGroup distribution:\n")
print(table(brca_prot_groups))

# -----------------------------------------------
# STEP 3: Handle NaN in normalization
# -----------------------------------------------
cat("\nSTEP 3: Normalizing with NaN handling...\n\n")

# Check for problematic values first
cat("Checking raw values...\n")
cat(sprintf("  Infinite values: %d\n",
            sum(is.infinite(
              as.matrix(brca_protein_raw)))))
cat(sprintf("  NaN values     : %d\n",
            sum(is.nan(
              as.matrix(brca_protein_raw)))))
cat(sprintf("  NA values      : %d\n",
            sum(is.na(brca_protein_raw))))
cat(sprintf("  Zero values    : %d\n",
            sum(brca_protein_raw == 0,
                na.rm=TRUE)))

# Replace infinite and NaN with NA
prot_mat <- as.matrix(brca_protein_raw)
prot_mat[is.infinite(prot_mat)] <- NA
prot_mat[is.nan(prot_mat)]      <- NA

# Remove proteins with >50% missing
miss_rate <- rowMeans(is.na(prot_mat))
prot_clean <- prot_mat[miss_rate < 0.50, ]

cat(sprintf("\nProteins after QC: %d\n",
            nrow(prot_clean)))

# Impute remaining missing
for(col in colnames(prot_clean)){
  idx <- is.na(prot_clean[, col])
  if(any(idx)){
    min_val <- min(prot_clean[, col],
                   na.rm=TRUE)
    prot_clean[idx, col] <- min_val / 2
  }
}

# Check if data needs log transformation
val_range <- range(prot_clean, na.rm=TRUE)
cat(sprintf("\nValue range: %.3f to %.3f\n",
            val_range[1], val_range[2]))

# LinkedOmics iTRAQ is already log-ratio
# No log transformation needed
# Just median center

if(abs(val_range[1]) < 20 &
   abs(val_range[2]) < 20){
  cat("Data appears to be log-ratio already\n")
  cat("Applying median centering only...\n\n")

  sample_medians    <- apply(prot_clean, 2,
                              median, na.rm=TRUE)
  brca_protein_norm <- sweep(prot_clean, 2,
                              sample_medians, "-")

} else {
  cat("Applying log2 transform + median center\n\n")

  # Shift to positive if needed
  if(val_range[1] <= 0){
    prot_clean <- prot_clean -
                  min(prot_clean, na.rm=TRUE) + 1
  }

  prot_log          <- log2(prot_clean)
  sample_medians    <- apply(prot_log, 2,
                              median, na.rm=TRUE)
  brca_protein_norm <- sweep(prot_log, 2,
                              sample_medians, "-")
}

# Final NaN check
nan_count <- sum(is.nan(brca_protein_norm))
if(nan_count > 0){
  brca_protein_norm[is.nan(
    brca_protein_norm)] <- NA
  cat(sprintf("Replaced %d NaN with NA\n\n",
              nan_count))
}

cat(sprintf("✅ Normalized!\n"))
cat(sprintf("  Proteins: %d\n",
            nrow(brca_protein_norm)))
cat(sprintf("  Samples : %d\n\n",
            ncol(brca_protein_norm)))

# Preview
cat("Preview (first 3 x 3):\n")
print(round(brca_protein_norm[1:3, 1:3], 4))

# -----------------------------------------------
# STEP 4: Create proper metadata
# -----------------------------------------------
cat("\nSTEP 4: Creating metadata...\n\n")

brca_prot_meta <- data.frame(
  sample_name = colnames(brca_protein_norm),
  group       = brca_prot_groups[
    colnames(brca_protein_raw) %in%
    colnames(brca_protein_norm)],
  stringsAsFactors = FALSE
)

cat("Group distribution in normalized data:\n")
print(table(brca_prot_meta$group))

# -----------------------------------------------
# STEP 5: Boxplot check
# -----------------------------------------------
cat("\nSTEP 5: Quality plot...\n\n")

library(ggplot2)

n_plot    <- min(20, ncol(brca_protein_norm))
plot_data <- data.frame(
  value  = as.vector(brca_protein_norm[, 1:n_plot]),
  sample = rep(colnames(brca_protein_norm)[1:n_plot],
               each=nrow(brca_protein_norm))
)

p <- ggplot(plot_data,
            aes(x=sample, y=value)) +
  geom_boxplot(fill="lightblue",
               outlier.size=0.5) +
  geom_hline(yintercept=0,
             color="red", linetype="dashed") +
  theme_bw() +
  theme(axis.text.x=element_text(
    angle=90, size=6)) +
  labs(title="Breast Proteomics: Normalized",
       y="Log Ratio Intensity")

print(p)
ggsave("04_Figures/Boxplot_Breast_Protein_Norm.pdf",
       p, width=10, height=6)
cat("✅ Boxplot saved\n")

# -----------------------------------------------
# SAVE
# -----------------------------------------------
cat("\nSaving...\n\n")

write.csv(
  brca_protein_norm,
  "02_ProcessedData/Breast/brca_protein_normalized.csv"
)
saveRDS(
  brca_prot_meta,
  "02_ProcessedData/Breast/brca_prot_meta.rds"
)
saveRDS(
  brca_prot_groups,
  "01_RawData/Breast/brca_prot_groups.rds"
)

size1 <- round(file.size(
  "02_ProcessedData/Breast/brca_protein_normalized.csv"
)/1024/1024, 2)

cat(sprintf(
  "✅ SAVED: brca_protein_normalized.csv  %.2f MB\n",
  size1))
cat("✅ SAVED: brca_prot_meta.rds\n")

# -----------------------------------------------
# SUMMARY
# -----------------------------------------------
cat("\n========================================\n")
cat("   STEP 2.3 SUMMARY\n")
cat("========================================\n\n")

cat(sprintf("Proteins         : %d\n",
            nrow(brca_protein_norm)))
cat(sprintf("Samples          : %d\n",
            ncol(brca_protein_norm)))
cat(sprintf("Cancer           : %d\n",
  sum(brca_prot_meta$group ==
      "EarlyBreastCancer")))
cat(sprintf("Normal           : %d\n",
  sum(brca_prot_meta$group == "Normal")))

missing_final <- round(
  100 * sum(is.na(brca_protein_norm)) /
  (nrow(brca_protein_norm) *
   ncol(brca_protein_norm)), 1)

cat(sprintf("Missing values   : %.1f%%\n",
            missing_final))
cat(sprintf("File size        : %.2f MB\n\n",
            size1))

cat("✅ STEP 2.3 COMPLETE\n")
cat("Continue to Step 2.4\n")
cat("========================================\n")