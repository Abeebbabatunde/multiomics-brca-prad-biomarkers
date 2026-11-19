# 01_stability_selection_rebuild.R
# Step 12.1 (Objectives 1 and 2): re-derive the panel for each layer by STABILITY selection inside patient-grouped resampling.
# The full selection (paired filter, LASSO, SVM-RFE, Boruta, consensus >= 2 of 3) is repeated on many random subsets of PATIENTS
# (a patient's tumour and normal samples always move together). The selection frequency of every feature across the subsets
# is the stability measure. Pre-specified rule: a feature enters the rebuilt panel if it is chosen in >= panel_freq of the subsets.
# Shuffled-label control: labels permuted WITHIN each patient; the same rule is applied, and the panel should be empty or tiny.
# Inputs come from 00_define_cohorts_and_inputs.R (same processing in both cancers; RNA transformed without the labels).
# Output: <layer>_selection_frequency_rebuild.csv, <layer>_rebuilt_panel.csv (+ shuffled control files).
# Run time: roughly 1-3 minutes per subset and layer (see the nested analysis, 40 fits in 0.5-2 hours); set n_subsets and
# layers_to_run to suit your computer. Edit project_dir only.

project_dir <- "C:/Users/USER/Documents/MyThesis_MultiOmics"
setwd(project_dir)
fn <- c("12_rebuild_patient_grouped/rebuild_functions.R", "rebuild_functions.R")   # in the project folder, or copied straight into it
fn <- fn[file.exists(fn)][1]
if (is.na(fn)) stop("rebuild_functions.R not found. Copy it into ", project_dir, " (or into a sub-folder called 12_rebuild_patient_grouped).")
source(fn)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db", update = FALSE, ask = FALSE)

## ------------------------------ settings ------------------------------
n_subsets      <- 50      # patient subsets per layer (more = more precise frequencies; the thesis nested analysis used 40 fits)
subset_frac    <- 0.75    # share of patients in each subset (3 of 4, as in the 4 outer folds of the nested analysis)
panel_freq     <- 0.60    # a feature joins the panel if selected in >= 60% of subsets (same threshold as the thesis stability selection)
run_shuffled   <- TRUE; n_subsets_shuffled <- 20
layers_to_run  <- c("PRAD_methylation", "BRCA_methylation", "BRCA_RNA", "PRAD_RNA")

layers <- list(
  PRAD_methylation = list(X = "12_rebuild_inputs/X_methylation_PRAD_rebuild.rds", meta = "12_rebuild_inputs/meta_PRAD_methylation_rebuild.csv", transpose = FALSE, ensembl = FALSE),
  BRCA_methylation = list(X = "12_rebuild_inputs/X_methylation_BRCA_rebuild.rds", meta = "12_rebuild_inputs/meta_BRCA_methylation_rebuild.csv", transpose = FALSE, ensembl = FALSE),
  BRCA_RNA         = list(X = "12_rebuild_inputs/X_rna_BRCA_rebuild.rds", meta = "12_rebuild_inputs/meta_BRCA_rna_rebuild.csv", transpose = FALSE, ensembl = TRUE),
  PRAD_RNA         = list(X = "12_rebuild_inputs/X_rna_PRAD_rebuild.rds", meta = "12_rebuild_inputs/meta_PRAD_rna_rebuild.csv", transpose = FALSE, ensembl = TRUE))

## ------------------------------ one layer ------------------------------
stability_run <- function(X, y, patient, n_sub, shuffle = FALSE) {
  ids <- unique(patient); sel_log <- vector("list", n_sub)
  for (b in seq_len(n_sub)) {
    yy <- y
    if (shuffle) for (id in ids) { i <- which(patient == id); yy[i] <- sample(y[i]) }   # permute WITHIN patient
    pick <- sample(ids, floor(length(ids) * subset_frac)); tr <- patient %in% pick
    sel_log[[b]] <- tryCatch(select_consensus(X[tr, , drop = FALSE], yy[tr], patient[tr]),
                             error = function(e) { message("subset ", b, " failed: ", conditionMessage(e)); character(0) })
    message(sprintf("%s subset %d/%d: %d features selected", ifelse(shuffle, "shuffled", "real"), b, n_sub, length(sel_log[[b]])))
  }
  sel_log
}
freq_table <- function(sel_log, ensembl) {
  n <- length(sel_log); f <- sort(table(unlist(sel_log)) / n, decreasing = TRUE)
  fd <- data.frame(feature = as.character(names(f)), selection_frequency = as.numeric(f), stringsAsFactors = FALSE)   # also works if nothing was selected
  if (nrow(fd) > 0 && ensembl) fd$symbol <- as.character(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
        keys = sub("\\.[0-9]+$", "", fd$feature), keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"))
  fd
}

run_layer <- function(name) {
  cfg <- layers[[name]]; cat("\n=========", name, "=========\n")
  set.seed(2026)
  X <- readRDS(cfg$X); if (cfg$transpose) X <- t(X)
  meta <- read.csv(cfg$meta, stringsAsFactors = FALSE); meta <- meta[match(rownames(X), meta$sample), ]
  stopifnot(!anyNA(meta$sample))
  y <- factor(meta$class, levels = c("Normal", "Tumour")); patient <- as.character(meta$patient)
  stopifnot(all(table(patient, y) > 0))
  X <- X[, apply(X, 2, function(v) !anyNA(v) && var(v) > 0), drop = FALSE]      # label-free clean-up
  cat("Matrix:", nrow(X), "samples x", ncol(X), "features;", length(unique(patient)), "patients\n")

  real <- stability_run(X, y, patient, n_subsets)
  fd <- freq_table(real, cfg$ensembl)
  write.csv(fd, paste0(name, "_selection_frequency_rebuild.csv"), row.names = FALSE)
  panel <- fd[fd$selection_frequency >= panel_freq, ]
  write.csv(panel, paste0(name, "_rebuilt_panel.csv"), row.names = FALSE)
  cat("\nTop 20 by selection frequency:\n"); print(head(fd, 20), row.names = FALSE)
  cat("Rebuilt panel (frequency >=", panel_freq, "):", nrow(panel), "features\n")

  if (run_shuffled) {
    shuf <- stability_run(X, y, patient, n_subsets_shuffled, shuffle = TRUE)
    fs <- freq_table(shuf, cfg$ensembl)
    write.csv(fs, paste0(name, "_selection_frequency_SHUFFLED.csv"), row.names = FALSE)
    cat("Shuffled-label control: features at or above the threshold =", sum(fs$selection_frequency >= panel_freq),
        "; highest frequency =", if (nrow(fs)) max(fs$selection_frequency) else 0, "\n")
  }
  invisible(list(freq = fd, panel = panel))
}
results <- lapply(layers_to_run, run_layer); names(results) <- layers_to_run
saveRDS(results, "rebuild_01_results.rds")
