## ============================================================
## Step 2.2 - Breast Methylation Differential Analysis
## Low-RAM (8GB) version - chunked processing
## ============================================================
## WHY THIS CRASHES ON 8GB MACHINES (read before running):
## The beta matrix itself (e.g. 450k probes x ~100 samples) is only
## ~350-400MB as a plain matrix. What actually kills an 8GB laptop is
## almost always one of these:
##   1. Loading the full minfi GenomicRatioSet / MethylSet object
##      (carries Beta + M + CN + full annotation slots = several GB)
##   2. Using champ.DMP() / dmpFinder() which loads the full
##      IlluminaHumanMethylation450kanno package (~300-500MB) plus
##      internal copies
##   3. A growing object in a loop, e.g. result <- rbind(result, x)
##      inside a for-loop over 450,000 probes - R recopies the whole
##      object every iteration (classic O(n^2) memory blowup)
##   4. Keeping the raw + normalized + filtered matrices all in
##      memory at once instead of removing each with rm() + gc()
##
## This script avoids all four: it only ever touches a lightweight
## beta MATRIX (not a minfi object), pre-filters low-information
## probes, processes probes in chunks, and explicitly frees memory
## between chunks.
## ============================================================

library(limma)
library(matrixStats)

## ---- 0. Settings ----
input_rds   <- "data/breast_methylation_normalized_beta.rds"  # <- adjust path
output_dir  <- "results/step2_2_methylation_diff/"
chunk_size  <- 40000      # probes per chunk; lower to 20000 if still tight
var_quantile_cutoff <- 0.25  # drop bottom 25% lowest-variance probes
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## ---- 1. Load ONLY the beta matrix, not a minfi object ----
## If what you currently have saved is a GenomicRatioSet/MethylSet,
## run this ONCE separately, save just the matrix, then restart R
## before running the rest of this script:
##   grset <- readRDS("data/your_minfi_object.rds")
##   beta  <- minfi::getBeta(grset)
##   saveRDS(beta, "data/breast_methylation_normalized_beta.rds")
##   rm(grset); gc()

beta <- readRDS(input_rds)
beta <- as.matrix(beta)   # ensure plain matrix, not data.frame (data.frames cost more RAM)
gc()

## sample_info must be a data.frame with rownames matching colnames(beta)
## and a "group" column: "Tumor" / "Normal"
sample_info <- readRDS("data/sample_info.rds")   # <- adjust path
sample_info <- sample_info[colnames(beta), , drop = FALSE]
group <- factor(sample_info$group, levels = c("Normal", "Tumor"))

cat("Loaded matrix:", nrow(beta), "probes x", ncol(beta), "samples\n")

## ---- 2. Pre-filter low-information probes (vectorized, cheap) ----
## rowVars() from matrixStats is C-vectorized - far less RAM/time than apply()
probe_var <- rowVars(beta, na.rm = TRUE)
keep_var  <- probe_var > quantile(probe_var, var_quantile_cutoff, na.rm = TRUE)

probe_na  <- rowSums(is.na(beta))
keep_na   <- probe_na == 0   # drop any probe with missing values (fast diff test needs complete data)

keep <- keep_var & keep_na
beta <- beta[keep, , drop = FALSE]
rm(probe_var, keep_var, probe_na, keep_na, keep); gc()

cat("After filtering:", nrow(beta), "probes remain\n")

## ---- 3. Build design matrix ONCE (reused for every chunk) ----
design <- model.matrix(~group)

## ---- 4. Chunked differential testing ----
probe_ids   <- rownames(beta)
n_probes    <- length(probe_ids)
chunk_starts <- seq(1, n_probes, by = chunk_size)

results_list <- vector("list", length(chunk_starts))

for (i in seq_along(chunk_starts)) {
  start_idx <- chunk_starts[i]
  end_idx   <- min(start_idx + chunk_size - 1, n_probes)

  chunk_mat <- beta[start_idx:end_idx, , drop = FALSE]

  fit  <- lmFit(chunk_mat, design)
  fit  <- eBayes(fit)
  tt   <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  tt$Probe <- rownames(chunk_mat)

  results_list[[i]] <- tt[, c("Probe", "logFC", "P.Value")]  # keep only what's needed

  ## checkpoint to disk in case the session dies mid-run
  saveRDS(results_list[1:i], file.path(output_dir, "checkpoint_partial_results.rds"))

  cat(sprintf("Chunk %d/%d done (probes %d-%d)\n", i, length(chunk_starts), start_idx, end_idx))

  rm(chunk_mat, fit, tt); gc()
}

## ---- 5. Combine + GLOBAL multiple-testing correction ----
## Important: do NOT just stack each chunk's adjusted p-values - BH/FDR
## correction must be recomputed across ALL probes together, not per chunk.
diff_results <- do.call(rbind, results_list)
diff_results$adj.P.Val <- p.adjust(diff_results$P.Value, method = "BH")
diff_results <- diff_results[order(diff_results$P.Value), ]

saveRDS(diff_results, file.path(output_dir, "breast_methylation_diffmeth_results.rds"))
write.csv(diff_results, file.path(output_dir, "breast_methylation_diffmeth_results.csv"), row.names = FALSE)

cat("Done. Significant probes (adj.P.Val < 0.05):", sum(diff_results$adj.P.Val < 0.05), "\n")

## ---- Extra tips if it still struggles on 8GB ----
## - Run via `Rscript breast_methylation_diff_lowRAM.R` from a terminal
##   instead of inside RStudio - RStudio itself holds 0.5-1GB+ overhead.
## - Close browser/Chrome tabs while running - each tab can use 200-500MB.
## - Drop chunk_size to 20000 or even 10000 if 40000 still hangs.
## - Annotate probes to genes AFTER this step, using a small probe->gene
##   lookup table (a few MB) rather than loading the full annotation
##   package again at this stage.