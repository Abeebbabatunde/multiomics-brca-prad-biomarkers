# ================================================
# STEP 3.2 — BREAST DMPs WITH CHECKPOINT SYSTEM
# Saves progress every 10,000 probes
# Can resume from where it stopped
# ================================================

cat("========================================\n")
cat("   STEP 3.2: BREAST DMPs\n")
cat("   WITH CHECKPOINT RECOVERY\n")
cat("========================================\n\n")

library(limma)

# -----------------------------------------------
# SETUP: Checkpoint folder
# -----------------------------------------------
checkpoint_dir <- "03_Results/Epigenomics/checkpoints"
dir.create(checkpoint_dir,
           recursive    = TRUE,
           showWarnings = FALSE)

checkpoint_file <- paste0(
  checkpoint_dir,
  "/brca_dmps_checkpoint.rds"
)
final_file <- paste0(
  checkpoint_dir,
  "/brca_dmps_all_results.rds"
)

cat(sprintf("Checkpoint folder: %s\n\n",
            checkpoint_dir))

# -----------------------------------------------
# STEP 1: Load methylation data
# -----------------------------------------------
cat("STEP 1: Loading methylation data...\n\n")

if(!exists("brca_meth_norm")){
  cat("Loading from file...\n")
  brca_meth_norm <- read.csv(
    "02_ProcessedData/Breast/brca_meth_normalized.csv",
    row.names   = 1,
    check.names = FALSE
  )
  cat(sprintf("✅ Loaded: %d probes x %d samples\n\n",
              nrow(brca_meth_norm),
              ncol(brca_meth_norm)))
}

if(!exists("meth_meta_brca")){
  meth_meta_brca <- readRDS(
    "02_ProcessedData/Breast/brca_meth_meta.rds")
}

cat(sprintf("Probes  : %d\n", nrow(brca_meth_norm)))
cat(sprintf("Samples : %d\n\n", ncol(brca_meth_norm)))

# -----------------------------------------------
# STEP 2: Prepare design matrix
# -----------------------------------------------
cat("STEP 2: Preparing design matrix...\n\n")

# Get group labels
group_meth <- factor(
  meth_meta_brca$group,
  levels = c("Normal","EarlyBreastCancer")
)

cat("Group distribution:\n")
print(table(group_meth))

design_meth <- model.matrix(~ 0 + group_meth)
colnames(design_meth) <- c("Normal",
                             "EarlyBreastCancer")

contrast_meth <- makeContrasts(
  EarlyBreastCancer - Normal,
  levels = design_meth
)

# -----------------------------------------------
# STEP 3: Run in chunks with checkpoints
# -----------------------------------------------
cat("\nSTEP 3: Running analysis in chunks...\n\n")

chunk_size  <- 10000  # Probes per chunk
total_probes <- nrow(brca_meth_norm)
n_chunks    <- ceiling(total_probes / chunk_size)

cat(sprintf("Total probes : %d\n", total_probes))
cat(sprintf("Chunk size   : %d\n", chunk_size))
cat(sprintf("Total chunks : %d\n\n", n_chunks))

# Check for existing checkpoint
start_chunk <- 1
all_results <- list()

if(file.exists(checkpoint_file)){
  cat("✅ Checkpoint found! Resuming...\n\n")
  checkpoint <- readRDS(checkpoint_file)
  start_chunk <- checkpoint$next_chunk
  all_results <- checkpoint$results
  cat(sprintf("Resuming from chunk %d of %d\n",
              start_chunk, n_chunks))
  cat(sprintf("Already completed: %d probes\n\n",
              (start_chunk-1) * chunk_size))
} else {
  cat("No checkpoint found — starting fresh\n\n")
}

# Run chunks
start_time <- Sys.time()

for(chunk_i in start_chunk:n_chunks){

  # Calculate probe indices for this chunk
  probe_start <- (chunk_i - 1) * chunk_size + 1
  probe_end   <- min(chunk_i * chunk_size,
                      total_probes)

  cat(sprintf(
    "Chunk %d/%d: probes %d to %d...",
    chunk_i, n_chunks,
    probe_start, probe_end))

  # Extract chunk
  meth_chunk <- as.matrix(
    brca_meth_norm[probe_start:probe_end, ])

  # Run limma on this chunk
  tryCatch({

    fit_chunk  <- lmFit(meth_chunk, design_meth)
    fit_chunk2 <- contrasts.fit(
      fit_chunk, contrast_meth)
    fit_chunk2 <- eBayes(fit_chunk2)

    # Get results
    results_chunk <- topTable(
      fit_chunk2,
      adjust = "BH",
      number = Inf,
      sort.by = "none"
    )

    all_results[[chunk_i]] <- results_chunk

    # Time estimate
    elapsed <- as.numeric(
      difftime(Sys.time(), start_time,
               units="mins"))
    chunks_done  <- chunk_i - start_chunk + 1
    chunks_left  <- n_chunks - chunk_i
    time_per_chunk <- elapsed / chunks_done
    time_left    <- round(
      chunks_left * time_per_chunk, 1)

    cat(sprintf(
      " ✅ done (%.1f min left)\n",
      time_left))

  }, error = function(e){
    cat(sprintf(" ❌ Error: %s\n",
                conditionMessage(e)))
  })

  # Save checkpoint every chunk
  checkpoint <- list(
    next_chunk = chunk_i + 1,
    results    = all_results,
    timestamp  = Sys.time()
  )
  saveRDS(checkpoint, checkpoint_file)
}

# -----------------------------------------------
# STEP 4: Combine all chunk results
# -----------------------------------------------
cat("\nSTEP 4: Combining all results...\n\n")

# Remove empty chunks
all_results <- all_results[
  !sapply(all_results, is.null)]

cat(sprintf("Chunks with results: %d\n",
            length(all_results)))

# Combine
brca_dmps_all <- do.call(rbind, all_results)

cat(sprintf("Total probes analyzed: %d\n\n",
            nrow(brca_dmps_all)))

# Filter significant DMPs
brca_dmps <- brca_dmps_all[
  !is.na(brca_dmps_all$adj.P.Val) &
  brca_dmps_all$adj.P.Val < 0.05 &
  abs(brca_dmps_all$logFC) > 0.2, ]

cat(sprintf("Significant DMPs: %d\n",
            nrow(brca_dmps)))
cat(sprintf("  Hypermethylated: %d\n",
            sum(brca_dmps$logFC > 0)))
cat(sprintf("  Hypomethylated : %d\n",
            sum(brca_dmps$logFC < 0)))

# -----------------------------------------------
# STEP 5: Save final results
# -----------------------------------------------
cat("\nSTEP 5: Saving results...\n\n")

write.csv(
  brca_dmps,
  "03_Results/Epigenomics/brca_DMPs.csv"
)
write.csv(
  brca_dmps_all,
  "03_Results/Epigenomics/brca_DMPs_all.csv"
)
saveRDS(
  brca_dmps_all,
  final_file
)

size1 <- round(file.size(
  "03_Results/Epigenomics/brca_DMPs.csv"
)/1024/1024, 2)

cat(sprintf(
  "✅ SAVED: brca_DMPs.csv  %.2f MB\n",
  size1))
cat("✅ SAVED: brca_DMPs_all.csv\n")

# Clean up checkpoint files
# (only after successful completion)
if(file.exists(checkpoint_file)){
  file.remove(checkpoint_file)
  cat("✅ Checkpoint files cleaned up\n")
}

# -----------------------------------------------
# FINAL SUMMARY
# -----------------------------------------------
total_time <- round(as.numeric(
  difftime(Sys.time(), start_time,
           units="mins")), 1)

cat("\n========================================\n")
cat("   STEP 3.2 COMPLETE\n")
cat("========================================\n\n")
cat(sprintf("Total time     : %.1f minutes\n",
            total_time))
cat(sprintf("Probes tested  : %d\n",
            nrow(brca_dmps_all)))
cat(sprintf("Significant    : %d\n",
            nrow(brca_dmps)))
cat(sprintf("Hypermethylated: %d\n",
            sum(brca_dmps$logFC > 0)))
cat(sprintf("Hypomethylated : %d\n",
            sum(brca_dmps$logFC < 0)))

cat("\n✅ Continue to Step 3.3\n")
cat("========================================\n")