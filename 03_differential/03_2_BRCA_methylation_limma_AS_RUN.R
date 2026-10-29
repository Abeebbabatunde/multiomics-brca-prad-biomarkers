## ============================================================
## Step 2.2 - Breast Methylation Differential Analysis (Tumor vs Normal)
## Low-RAM (8GB) version - chunked processing
## ============================================================
## CONFIRMED STRUCTURE (from diagnostics):
##  - 01_RawData/Breast/brca_meth.rds is a RangedSummarizedExperiment
##    (TCGAbiolinks output), 485,577 probes x 649 samples, ~2.5GB as
##    a raw matrix, peaking ~5.3GB during extraction - this alone is
##    the main memory bottleneck on an 8GB machine.
##  - The unnamed assay ('') is the beta value matrix - extract with
##    assay(beta), not minfi::getBeta().
##  - rowData(beta)$MASK_general flags unreliable/cross-reactive probes
##    (Zhou et al. masking) - legitimate quality filter, not arbitrary.
##  - 01_RawData/Breast/brca_meth_meta.rds only has Sample/PatientID/
##    Group columns, and Group is uniform ("EarlyBreastCancer") for
##    all 649 samples - it does NOT distinguish tumor vs normal.
##  - Tumor/Normal status must be derived from the TCGA barcode itself
##    (segment 4, sample type code: 01 = Tumor, 11 = Normal).
## ============================================================

library(SummarizedExperiment)
library(limma)
library(matrixStats)

## ---- 0. Settings ----
input_rds  <- "01_RawData/Breast/brca_meth.rds"
checkpoint_dir <- "02_ProcessedData/Breast/"
output_dir <- "results/step2_2_methylation_diff/"
chunk_size <- 40000
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

## ---- 1. Load, extract matrix + mask flag, then immediately drop the
##         heavy S4 object ----
beta_obj     <- readRDS(input_rds)
mask_general <- rowData(beta_obj)$MASK_general
beta_matrix  <- assay(beta_obj)
rm(beta_obj); gc()

cat("Raw matrix:", nrow(beta_matrix), "probes x", ncol(beta_matrix), "samples\n")

## ---- 2. Derive Tumor/Normal group from the TCGA barcode ----
barcode_parts     <- strsplit(colnames(beta_matrix), "-")
sample_type_field  <- sapply(barcode_parts, function(x) x[4])
sample_type_code   <- substr(sample_type_field, 1, 2)

group <- factor(ifelse(sample_type_code == "11", "Normal",
                 ifelse(sample_type_code == "01", "Tumor", NA)),
                 levels = c("Normal", "Tumor"))

keep_samples <- !is.na(group)
beta_matrix  <- beta_matrix[, keep_samples]
group        <- group[keep_samples]

cat("After dropping non-Tumor/Normal barcodes:", ncol(beta_matrix), "samples\n")
print(table(group))

## ---- 3. Filter probes: drop masked (unreliable) + any-NA probes ----
keep_probes <- !mask_general & rowSums(is.na(beta_matrix)) == 0
beta_matrix <- beta_matrix[keep_probes, ]
rm(mask_general, keep_probes); gc()

cat("After MASK_general + NA filtering:", nrow(beta_matrix), "probes remain\n")

## ---- 4. Checkpoint the cleaned inputs before the long loop ----
## So a crash during the differential step doesn't force redoing
## the extraction + filtering above.
saveRDS(beta_matrix, file.path(checkpoint_dir, "brca_meth_filtered_matrix.rds"))
saveRDS(group, file.path(checkpoint_dir, "brca_meth_filtered_group.rds"))

## ---- 5. Build design matrix ONCE ----
design <- model.matrix(~group)

## ---- 6. Chunked differential testing ----
probe_ids    <- rownames(beta_matrix)
n_probes     <- length(probe_ids)
chunk_starts <- seq(1, n_probes, by = chunk_size)
results_list <- vector("list", length(chunk_starts))

for (i in seq_along(chunk_starts)) {
  start_idx <- chunk_starts[i]
  end_idx   <- min(start_idx + chunk_size - 1, n_probes)

  chunk_mat <- beta_matrix[start_idx:end_idx, , drop = FALSE]

  fit <- lmFit(chunk_mat, design)
  fit <- eBayes(fit)
  tt  <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  tt$Probe <- rownames(chunk_mat)

  results_list[[i]] <- tt[, c("Probe", "logFC", "P.Value")]

  saveRDS(results_list[1:i], file.path(output_dir, "checkpoint_partial_results.rds"))
  cat(sprintf("Chunk %d/%d done (probes %d-%d)\n", i, length(chunk_starts), start_idx, end_idx))

  rm(chunk_mat, fit, tt); gc()
}

## ---- 7. Combine + GLOBAL multiple-testing correction ----
diff_results <- do.call(rbind, results_list)
diff_results$adj.P.Val <- p.adjust(diff_results$P.Value, method = "BH")
diff_results <- diff_results[order(diff_results$P.Value), ]

saveRDS(diff_results, file.path(output_dir, "breast_methylation_diffmeth_results.rds"))
write.csv(diff_results, file.path(output_dir, "breast_methylation_diffmeth_results.csv"), row.names = FALSE)

cat("Done. Significant probes (adj.P.Val < 0.05):", sum(diff_results$adj.P.Val < 0.05), "\n")