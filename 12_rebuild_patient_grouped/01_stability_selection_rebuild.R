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
## ---- selection functions (built into this file, so no separate functions file is needed) ----
# rebuild_functions.R
# Shared functions for the rebuild (Step 12). The selection functions are the same as in
# 11_pre_submission_reanalyses/01_nested_patient_cv.R, so the rebuild uses the thesis methods unchanged.
# Every function that uses the class labels works on ONE training set and is told the patient of each sample,
# so a patient's tumour and normal samples are never split.

suppressPackageStartupMessages({library(glmnet); library(limma); library(Boruta); library(e1071); library(pROC)})

patient_folds <- function(patient, k) {
  ids <- unique(patient); k <- min(k, length(ids))
  f <- sample(rep(seq_len(k), length.out = length(ids))); f[match(patient, ids)]
}
filter_features <- function(Xtr, ytr, ptr, n_keep) {          # paired moderated t, patient as blocking factor
  design <- model.matrix(~ factor(ptr) + ytr)
  fit <- eBayes(lmFit(t(Xtr), design))
  tt <- topTable(fit, coef = ncol(design), number = Inf, sort.by = "P")
  rownames(tt)[seq_len(min(n_keep, nrow(tt)))]
}
sel_lasso <- function(Xtr, ytr, ptr, B, freq) {               # LASSO + stability selection (subsample patients)
  ids <- unique(ptr)
  cv <- cv.glmnet(Xtr, ytr, family = "binomial", alpha = 1, foldid = patient_folds(ptr, min(5, length(ids))))
  lam <- cv$lambda.1se; hits <- setNames(numeric(ncol(Xtr)), colnames(Xtr))
  for (b in seq_len(B)) {
    sub <- sample(ids, floor(length(ids) / 2)); i <- ptr %in% sub
    if (length(unique(ytr[i])) < 2) next
    fit <- glmnet(Xtr[i, , drop = FALSE], ytr[i], family = "binomial", alpha = 1, lambda = lam)
    nz <- as.vector(coef(fit))[-1] != 0; hits[nz] <- hits[nz] + 1
  }
  names(hits)[hits / B >= freq]
}
svm_rfe <- function(Xs, ys, n_keep, drop = 0.2) {
  feats <- colnames(Xs)
  while (length(feats) > n_keep) {
    fit <- svm(Xs[, feats, drop = FALSE], ys, kernel = "linear", cost = 1, scale = TRUE)
    w <- abs(as.vector(t(fit$coefs) %*% fit$SV))
    n_drop <- min(max(1, floor(length(feats) * drop)), length(feats) - n_keep)
    feats <- feats[order(w, decreasing = TRUE)][seq_len(length(feats) - n_drop)]
  }
  feats
}
sel_svm <- function(Xtr, ytr, ptr, B = 20, freq, n_keep) {    # SVM-RFE + stability selection
  ids <- unique(ptr); hits <- setNames(numeric(ncol(Xtr)), colnames(Xtr))
  for (b in seq_len(B)) {
    sub <- sample(ids, floor(length(ids) * 0.7)); i <- ptr %in% sub
    if (length(unique(ytr[i])) < 2) next
    s <- svm_rfe(Xtr[i, , drop = FALSE], ytr[i], n_keep); hits[s] <- hits[s] + 1
  }
  names(hits)[hits / B >= freq]
}
sel_boruta <- function(Xtr, ytr) {
  res <- Boruta(as.data.frame(Xtr), ytr, maxRuns = 100, doTrace = 0)
  getSelectedAttributes(res, withTentative = FALSE)
}
# One full selection on one set of patients: filter -> LASSO, SVM-RFE, Boruta -> consensus (>= 2 of 3 methods)
sel_lasso_single <- function(Xtr, ytr, ptr) {                   # one LASSO fit at lambda.1se, no inner stability
  cv <- cv.glmnet(Xtr, ytr, family = "binomial", alpha = 1, foldid = patient_folds(ptr, min(5, length(unique(ptr)))))
  nz <- as.vector(coef(cv, s = "lambda.1se"))[-1] != 0
  colnames(Xtr)[nz]
}
# inner_stability = TRUE : thesis version (LASSO and SVM-RFE each use their own inner stability threshold)
# inner_stability = FALSE: one fit per selector; the stability is then measured ONCE, by the selection frequency
#                          across the outer patient subsets (avoids applying a stability threshold twice)
select_consensus <- function(Xtr, ytr, ptr, n_filter = 500, n_svm_keep = 20, stab_B = 50, stab_freq = 0.60, inner_stability = TRUE) {
  keep <- filter_features(Xtr, ytr, ptr, n_filter); Xf <- Xtr[, keep, drop = FALSE]
  if (inner_stability) {
    s_l <- sel_lasso(Xf, ytr, ptr, stab_B, stab_freq)
    s_s <- sel_svm(Xf, ytr, ptr, B = 20, freq = stab_freq, n_keep = n_svm_keep)
  } else {
    s_l <- sel_lasso_single(Xf, ytr, ptr)
    s_s <- svm_rfe(Xf, ytr, n_svm_keep)
  }
  s_b <- sel_boruta(Xf, ytr)
  votes <- table(factor(c(s_l, s_s, s_b), levels = keep))
  names(votes)[votes >= 2]
}

cat("Script version: self-contained, inner_stability option present\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db", update = FALSE, ask = FALSE)

## ------------------------------ settings ------------------------------
n_subsets      <- 50      # patient subsets per layer (more = more precise frequencies; the thesis nested analysis used 40 fits)
subset_frac    <- 0.75    # share of patients in each subset (3 of 4, as in the 4 outer folds of the nested analysis)
panel_freq     <- 0.60    # a feature joins the panel if selected in >= 60% of subsets (same threshold as the thesis stability selection)
run_shuffled   <- TRUE; n_subsets_shuffled <- 20
inner_stability <- FALSE   # FALSE = stability measured once (outer subsets); TRUE = thesis version (stability inside each selector as well)
layers_to_run  <- c("PRAD_methylation", "BRCA_methylation", "BRCA_RNA", "PRAD_RNA")   # FULL RUN (pre-specified: 50 subsets, 20 shuffled, 0.60, 75% of patients)

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
    sel_log[[b]] <- tryCatch(select_consensus(X[tr, , drop = FALSE], yy[tr], patient[tr], inner_stability = inner_stability),
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
  cfg <- layers[[name]]; cat("\n=========", name, "=========\n"); t_layer <- Sys.time()
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
  cat(sprintf("Time for %s: %.1f minutes\n", name, as.numeric(difftime(Sys.time(), t_layer, units = "mins"))))
  invisible(list(freq = fd, panel = panel))
}
results <- lapply(layers_to_run, run_layer); names(results) <- layers_to_run
saveRDS(results, "rebuild_01_results.rds")
