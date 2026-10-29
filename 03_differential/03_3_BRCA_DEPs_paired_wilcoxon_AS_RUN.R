# ============================================================
> 
> dep_results$direction <- ifelse(dep_results$mean_diff > 0, "Up in Tumor", "Down in Tumor")
> 
> cat("\nSTEP 4: Significant DEPs found\n")# ============================================================
# STEP 3.3: Breast DEPs (Differential Expression Proteins)
# Paired tumor vs normal analysis on normalized RPPA data
# ============================================================

library(dplyr)
library(tidyr)

# ---- 1. Load Step 2.3 output ----
brca_protein_normalized <- read.csv("02_Normalized/Breast/brca_protein_normalized.csv", 
                                      stringsAsFactors = FALSE)

# ---- 2. Pivot to wide format: one row per protein, columns = patients, split by tumor/normal ----
# Add patient_id column (strip sample suffix to get the person)
brca_protein_normalized$patient_id <- substr(brca_protein_normalized$sample_barcode, 1, 12)

# Separate tumor and normal values, keyed by protein + patient
tumor_wide <- brca_protein_normalized %>%
  filter(sample_type == "Primary Tumor") %>%
  select(peptide_target, patient_id, value) %>%
  rename(tumor_value = value)

normal_wide <- brca_protein_normalized %>%
  filter(sample_type == "Solid Tissue Normal") %>%
  select(peptide_target, patient_id, value) %>%
  rename(normal_value = value)

# Join tumor and normal values for the same patient + protein (this enforces pairing)
paired_data <- inner_join(tumor_wide, normal_wide, by = c("peptide_target", "patient_id"))

cat("Paired data check — should be 457 proteins × 30 patients =", 457*30, "rows\n")
cat("Actual rows:", nrow(paired_data), "\n")

# ---- 3. Run paired Wilcoxon signed-rank test per protein ----
# (Wilcoxon chosen over paired t-test since RPPA log-ratios can have outliers/non-normality;
#  switch to t.test(paired=TRUE) below if you confirm normality per-protein instead)

cat("\nSTEP 3: Running paired differential tests per protein\n")

dep_results <- paired_data %>%
  group_by(peptide_target) %>%
  summarise(
    mean_tumor   = mean(tumor_value, na.rm = TRUE),
    mean_normal  = mean(normal_value, na.rm = TRUE),
    mean_diff    = mean(tumor_value - normal_value, na.rm = TRUE),
    p_value      = tryCatch(
      wilcox.test(tumor_value, normal_value, paired = TRUE)$p.value,
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  arrange(p_value)

# ---- 4. Multiple testing correction (FDR / Benjamini-Hochberg) ----
dep_results$p_adj <- p.adjust(dep_results$p_value, method = "fdr")

# ---- 5. Flag significant DEPs ----
# Common thresholds: p_adj < 0.05, and a minimum effect size (|mean_diff| > 0.5 is a common RPPA cutoff,
# but adjust based on your field's convention)
sig_threshold_padj <- 0.05
sig_threshold_diff <- 0.5

dep_results$significant <- dep_results$p_adj < sig_threshold_padj & 
                            abs(dep_results$mean_diff) > sig_threshold_diff

dep_results$direction <- ifelse(dep_results$mean_diff > 0, "Up in Tumor", "Down in Tumor")

cat("\nSTEP 4: Significant DEPs found\n")
cat("Total proteins tested:", nrow(dep_results), "\n")
cat("Significant DEPs (p_adj <", sig_threshold_padj, "and |diff| >", sig_threshold_diff, "):", 
    sum(dep_results$significant), "\n")

print(head(dep_results[dep_results$significant, ], 20))

# ---- 6. Save Step 3.3 output ----
deg_dir <- "03_DifferentialExpression/Breast"  # adjust to match your actual pipeline folder naming
dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(dep_results, file.path(deg_dir, "brca_DEPs.csv"), row.names = FALSE)

cat("\n✅ Saved:", file.path(deg_dir, "brca_DEPs.csv"), "\n")
cat("=== STEP 3.3 COMPLETE ===\n")
