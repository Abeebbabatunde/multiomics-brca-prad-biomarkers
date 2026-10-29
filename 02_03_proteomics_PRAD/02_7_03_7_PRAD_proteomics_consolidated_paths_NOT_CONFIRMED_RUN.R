## =========================================================================
## MyThesis_MultiOmics — Prostate Proteomics
## Step 2.7: Protein QC / Filtering, Log2 Transformation, Normalization
## Step 3.7: Differential Protein Expression (Localized Cancer vs Adjacent Normal)
##
## Data source : E-PROT-21
##   28 Localized prostate cancer
##    8 Adjacent non-malignant
## Input file  : prad_protein_raw.tsv
## =========================================================================

## ---- 0. Packages --------------------------------------------------------
required_pkgs <- c("limma", "edgeR", "preprocessCore", "dplyr", "tibble",
                    "ggplot2", "ggrepel", "pheatmap")

new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(new_pkgs, update = FALSE, ask = FALSE)
}

library(limma)
library(edgeR)
library(preprocessCore)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(pheatmap)

set.seed(2026)

## ---- 1. Paths -------------------------------------------------------------
data_dir   <- "01_RawData/Prostate"
out_dir_processed <- "02_ProcessedData/Prostate"   # Step 2.7 outputs (normalized data)
out_dir_results   <- "03_Results/Proteomics"       # Step 3.7 outputs (DEP results)
dir.create(out_dir_processed, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_results, recursive = TRUE, showWarnings = FALSE)

raw_file <- file.path(data_dir, "prad_protein_raw.tsv")

## =========================================================================
## STEP 2.7 — Prostate Proteomics: QC / Filtering, Log2 Transform, Normalization
## =========================================================================

## ---- 2.7.1 Load raw protein abundance matrix ------------------------------
# This is an Expression Atlas export (E-PROT-21). The file has 4 '#' comment
# lines before the real header row, and TWO identifier columns:
#   Gene ID   (Ensembl gene ID, unique -> used as rownames)
#   Gene Name (HGNC symbol, kept as an annotation lookup for later reporting)
# Sample columns are labeled e.g. "individual 1, neoplasm" or
# "individual 1, normal tissue adjacent to tumor" -- no name-guessing needed.

raw <- read.delim(raw_file, header = TRUE, sep = "\t", skip = 4,
                   check.names = FALSE, stringsAsFactors = FALSE,
                   comment.char = "")

stopifnot(all(c("Gene ID", "Gene Name") %in% colnames(raw)))

# Keep Gene ID <-> Gene Name mapping for later annotation of DEP results
gene_annotation <- raw %>% dplyr::select(`Gene ID`, `Gene Name`)

prot_mat <- raw %>%
  dplyr::select(-`Gene Name`) %>%
  column_to_rownames(var = "Gene ID") %>%
  as.matrix()

storage.mode(prot_mat) <- "numeric"

cat("Loaded matrix:", nrow(prot_mat), "proteins x", ncol(prot_mat), "samples\n")

## ---- 2.7.2 Build sample metadata (E-PROT-21 design) -----------------------
# Group is read directly off the column labels -- Expression Atlas already
# encodes it explicitly ("neoplasm" vs "normal tissue adjacent to tumor"),
# so no pattern-guessing or external metadata file is needed here.

sample_names <- colnames(prot_mat)

group <- case_when(
  grepl("normal tissue adjacent to tumor", sample_names, ignore.case = TRUE) ~ "Adjacent_Normal",
  grepl("neoplasm", sample_names, ignore.case = TRUE) ~ "Localized_Cancer",
  TRUE ~ NA_character_
)

if (any(is.na(group))) {
  stop("Some sample columns didn't match 'neoplasm' or 'normal tissue adjacent to tumor'. ",
       "Unmatched columns: ", paste(sample_names[is.na(group)], collapse = "; "))
}

meta <- data.frame(sample = sample_names, group = group, stringsAsFactors = FALSE)
meta$group <- factor(meta$group, levels = c("Adjacent_Normal", "Localized_Cancer"))

stopifnot(sum(meta$group == "Localized_Cancer")   == 28)
stopifnot(sum(meta$group == "Adjacent_Normal")    == 8)

table(meta$group)

## ---- 2.7.3 Protein-level QC / filtering -----------------------------------
# a) Drop proteins with excessive missingness (keep proteins with valid
#    (non-NA, non-zero) values in >= 70% of samples)
missing_thresh <- 0.70
keep_valid <- rowMeans(!is.na(prot_mat) & prot_mat != 0) >= missing_thresh
prot_filt <- prot_mat[keep_valid, ]

cat(sum(keep_valid), "of", nrow(prot_mat),
    "proteins retained after missingness filter (>=", missing_thresh * 100, "% valid)\n")

# b) Drop near-zero-variance proteins (uninformative across samples)
row_var <- apply(prot_filt, 1, var, na.rm = TRUE)
keep_var <- row_var > 0 & !is.na(row_var)
prot_filt <- prot_filt[keep_var, ]

cat(nrow(prot_filt), "proteins retained after variance filter\n")

## ---- 2.7.4 Impute remaining missing values --------------------------------
# Simple approach: replace remaining NAs with per-protein minimum
# observed value (common for MNAR-type missingness in MS proteomics).
# Swap for a more rigorous method (e.g. QRILC, KNN) if required.
impute_min <- function(x) {
  if (any(is.na(x))) x[is.na(x)] <- min(x, na.rm = TRUE) * 0.5
  x
}
prot_filt <- t(apply(prot_filt, 1, impute_min))

## ---- 2.7.5 Log2 transformation --------------------------------------------
# Only log-transform if data is not already logged (heuristic: max value > 50)
if (max(prot_filt, na.rm = TRUE) > 50) {
  prot_log2 <- log2(prot_filt + 1)
} else {
  prot_log2 <- prot_filt
  message("Data appears already log-scaled; skipping additional log2 transform.")
}

## ---- 2.7.6 Normalization ---------------------------------------------------
# Quantile normalization across samples (standard for MS/RPPA-style
# proteomics matrices; makes sample distributions comparable)
prot_norm <- normalize.quantiles(prot_log2)
dimnames(prot_norm) <- dimnames(prot_log2)

## ---- 2.7.7 QC plots --------------------------------------------------------
pdf(file.path(out_dir_processed, "prad_protein_QC_boxplots.pdf"), width = 10, height = 6)
par(mfrow = c(1, 2))
boxplot(prot_log2, las = 2, main = "Before normalization", ylab = "log2 abundance",
        col = ifelse(meta$group == "Localized_Cancer", "tomato", "steelblue"))
boxplot(prot_norm, las = 2, main = "After quantile normalization", ylab = "log2 abundance",
        col = ifelse(meta$group == "Localized_Cancer", "tomato", "steelblue"))
dev.off()

## PCA sanity check
pca <- prcomp(t(prot_norm), scale. = TRUE)
pca_df <- data.frame(pca$x[, 1:2], group = meta$group, sample = meta$sample)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = group, label = sample)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PRAD Proteomics — PCA (post-normalization)",
       x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)"))
ggsave(file.path(out_dir_processed, "prad_protein_PCA.pdf"), p_pca, width = 7, height = 5)

## ---- 2.7.8 Save normalized matrix ------------------------------------------
# Primary output — matches the exact filename checked by the progress tracker
# (Step 2.7, "Prostate protein normalized" -> 02_ProcessedData/Prostate/prad_protein_normalized.csv)
write.csv(prot_norm, file.path(out_dir_processed, "prad_protein_normalized.csv"))

# Supplementary outputs
saveRDS(prot_norm, file.path(out_dir_processed, "prad_protein_normalized.rds"))
write.table(prot_norm, file.path(out_dir_processed, "prad_protein_normalized.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
saveRDS(meta, file.path(out_dir_processed, "prad_protein_metadata.rds"))

cat("\n[Step 2.7 complete] Normalized matrix:", nrow(prot_norm), "proteins x",
    ncol(prot_norm), "samples\n")

## =========================================================================
## STEP 3.7 — Prostate DEPs: Differential Protein Expression
## Localized Cancer vs Adjacent Normal
## =========================================================================

## ---- 3.7.1 Design matrix ---------------------------------------------------
design <- model.matrix(~ 0 + group, data = meta)
colnames(design) <- levels(meta$group)   # Adjacent_Normal, Localized_Cancer

## ---- 3.7.2 Fit limma model --------------------------------------------------
fit <- lmFit(prot_norm, design)

contrast_matrix <- makeContrasts(
  Cancer_vs_Normal = Localized_Cancer - Adjacent_Normal,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, trend = TRUE)   # trend=TRUE recommended for proteomics/RPPA-style data

## ---- 3.7.3 Extract results --------------------------------------------------
deps_all <- topTable(fit2, coef = "Cancer_vs_Normal", number = Inf, sort.by = "P") %>%
  rownames_to_column(var = "Protein") %>%
  rename(log2FC = logFC, p_value = P.Value, FDR = adj.P.Val) %>%
  left_join(gene_annotation, by = c("Protein" = "Gene ID")) %>%
  rename(GeneName = `Gene Name`) %>%
  relocate(GeneName, .after = Protein)

## ---- 3.7.4 Define significance thresholds -----------------------------------
logfc_cutoff <- 1        # |log2FC| >= 1 (2-fold)
fdr_cutoff   <- 0.05

deps_all <- deps_all %>%
  mutate(
    Significance = case_when(
      FDR < fdr_cutoff & log2FC >=  logfc_cutoff ~ "Up in Cancer",
      FDR < fdr_cutoff & log2FC <= -logfc_cutoff ~ "Down in Cancer",
      TRUE ~ "Not Significant"
    )
  )

deps_sig <- deps_all %>% filter(Significance != "Not Significant")

cat("\n[Step 3.7] Differential protein expression (Localized Cancer vs Adjacent Normal)\n")
cat("Total proteins tested :", nrow(deps_all), "\n")
cat("Significant DEPs      :", nrow(deps_sig),
    "(", sum(deps_sig$Significance == "Up in Cancer"), "up /",
    sum(deps_sig$Significance == "Down in Cancer"), "down )\n")

## ---- 3.7.5 Save DEP tables --------------------------------------------------
# Primary output — matches the exact filename checked by the progress tracker
# (Step 3.7, "Prostate DEPs" -> 03_Results/Proteomics/prad_DEPs.csv)
write.csv(deps_all, file.path(out_dir_results, "prad_DEPs.csv"), row.names = FALSE)

# Supplementary outputs (kept for convenience, not required by the tracker)
write.table(deps_all, file.path(out_dir_results, "prad_DEPs_all.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(deps_sig, file.path(out_dir_results, "prad_DEPs_significant.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(deps_all, file.path(out_dir_results, "prad_DEPs_all.rds"))

## ---- 3.7.6 Volcano plot -----------------------------------------------------
p_volcano <- ggplot(deps_all, aes(x = log2FC, y = -log10(p_value), color = Significance)) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_color_manual(values = c(
    "Up in Cancer"     = "firebrick",
    "Down in Cancer"   = "steelblue",
    "Not Significant"  = "grey70"
  )) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed", color = "grey40") +
  geom_text_repel(
    data = deps_sig %>% arrange(FDR) %>% head(15),
    aes(label = Protein), size = 3, max.overlaps = 20, show.legend = FALSE
  ) +
  theme_minimal(base_size = 12) +
  labs(title = "PRAD DEPs — Localized Cancer vs Adjacent Normal",
       subtitle = paste0("|log2FC| >= ", logfc_cutoff, ", FDR < ", fdr_cutoff),
       x = "log2 Fold Change", y = "-log10(p-value)")

ggsave(file.path(out_dir_results, "prad_DEP_volcano_plot.pdf"), p_volcano, width = 7, height = 6)

## ---- 3.7.7 Heatmap of top DEPs ----------------------------------------------
top_n <- min(50, nrow(deps_sig))

if (top_n > 1) {
  top_proteins <- deps_sig %>% arrange(FDR) %>% slice_head(n = top_n) %>% pull(Protein)

  heat_mat <- prot_norm[top_proteins, , drop = FALSE]
  heat_mat_z <- t(scale(t(heat_mat)))   # row z-score for visualization

  annotation_col <- data.frame(Group = meta$group)
  rownames(annotation_col) <- meta$sample

  pheatmap(
    heat_mat_z,
    annotation_col = annotation_col,
    show_rownames = ifelse(top_n <= 50, TRUE, FALSE),
    show_colnames = TRUE,
    fontsize_row = 6,
    color = colorRampPalette(c("steelblue", "white", "firebrick"))(100),
    main = paste0("Top ", top_n, " Prostate DEPs (Cancer vs Adjacent Normal)"),
    filename = file.path(out_dir_results, "prad_DEP_heatmap.pdf"),
    width = 8, height = 10
  )
} else {
  message("Fewer than 2 significant DEPs — skipping heatmap.")
}

cat("\n[Step 3.7 complete] Outputs written to:", out_dir_results, "\n")