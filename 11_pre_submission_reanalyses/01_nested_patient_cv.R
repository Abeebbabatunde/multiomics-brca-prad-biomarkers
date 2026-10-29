# 01_nested_patient_cv.R
# FIX 1: nested, patient-grouped resampling for the four layers (breast and prostate; methylation and RNA).
# Every step that uses the class labels (paired moderated-t filter, LASSO, SVM-RFE and Boruta selection,
# consensus, final model) is repeated inside each training fold, and folds are formed by PATIENT.
# A shuffled-label control permutes the labels WITHIN each patient.
# Requires the files made by 00_build_inputs.R.  Run time: about 0.5-2 hours per layer.
# NOTE: stats::predict is called explicitly because mixOmics (used in script 03a) masks predict().

project_dir <- "C:/Users/USER/Documents/MyThesis_MultiOmics"
setwd(project_dir)
for (p in c("glmnet", "Boruta", "e1071", "pROC")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in c("limma", "org.Hs.eg.db")) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)
suppressPackageStartupMessages({library(glmnet); library(limma); library(Boruta); library(e1071); library(pROC)})

## ------------------------------ settings ------------------------------
n_outer <- 4; n_reps <- 10; n_reps_shuffle <- 5
n_filter <- 500; n_svm_keep <- 20; stab_B <- 50; stab_freq <- 0.60
layers_to_run <- c("PRAD_methylation", "BRCA_methylation", "BRCA_RNA", "PRAD_RNA")

layers <- list(
  PRAD_methylation = list(X = "X_methylation_PRAD.rds", meta = "meta_PRAD.csv",
                          prefix = "nested_PRAD_methylation_FINAL",
                          panel_file = "10_FinalResults/Step10_FinalPanel/Step10_2_PRAD_final_panel.csv", omics = "methylation", ensembl = FALSE),
  BRCA_methylation = list(X = "X_methylation_BRCA.rds", meta = "meta_BRCA.csv",
                          prefix = "nested_BRCA_methylation_FINAL",
                          panel_file = "10_FinalResults/Step10_FinalPanel/Step10_1_BRCA_final_panel.csv", omics = "methylation", ensembl = FALSE),
  BRCA_RNA         = list(X = "X_rna_BRCA.rds", meta = "meta_BRCA_rna.csv",
                          prefix = "nested_BRCA_RNA_FINAL",
                          panel_file = "10_FinalResults/Step10_FinalPanel/Step10_1_BRCA_final_panel.csv", omics = "RNA", ensembl = TRUE),
  PRAD_RNA         = list(X = "PRAD_rna_log2.rds", meta = "meta_PRAD_matched.csv",
                          prefix = "nested_PRAD_RNA_log2", transpose = TRUE,   # PRAD_rna_log2.rds is genes x samples
                          panel_file = "10_FinalResults/Step10_FinalPanel/Step10_2_PRAD_final_panel.csv", omics = "RNA", ensembl = TRUE))

## ------------------------------ functions ------------------------------
patient_folds <- function(patient, k) {
  ids <- unique(patient); k <- min(k, length(ids))
  f <- sample(rep(seq_len(k), length.out = length(ids))); f[match(patient, ids)]
}
filter_features <- function(Xtr, ytr, ptr, n_keep) {                       # paired moderated t, patient as blocking factor
  design <- model.matrix(~ factor(ptr) + ytr)
  fit <- eBayes(lmFit(t(Xtr), design))
  tt <- topTable(fit, coef = ncol(design), number = Inf, sort.by = "P")
  rownames(tt)[seq_len(min(n_keep, nrow(tt)))]
}
sel_lasso <- function(Xtr, ytr, ptr, B, freq) {                            # LASSO + stability selection (subsample patients)
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
sel_svm <- function(Xtr, ytr, ptr, B = 20, freq, n_keep) {                 # SVM-RFE + stability selection
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
fit_predict <- function(Xtr, ytr, ptr, Xte) {                              # everything inside ONE training set
  keep <- filter_features(Xtr, ytr, ptr, n_filter); Xf <- Xtr[, keep, drop = FALSE]
  s_l <- sel_lasso(Xf, ytr, ptr, stab_B, stab_freq)
  s_s <- sel_svm(Xf, ytr, ptr, B = 20, freq = stab_freq, n_keep = n_svm_keep)
  s_b <- sel_boruta(Xf, ytr)
  votes <- table(factor(c(s_l, s_s, s_b), levels = keep)); sel <- names(votes)[votes >= 2]   # consensus: >= 2 of 3 methods
  if (length(sel) < 2) sel <- if (length(s_l) >= 2) s_l else keep[1:5]
  cv <- cv.glmnet(Xtr[, sel, drop = FALSE], ytr, family = "binomial", alpha = 0,             # ridge final classifier
                  foldid = patient_folds(ptr, min(5, length(unique(ptr)))))
  prob <- as.numeric(stats::predict(cv, Xte[, sel, drop = FALSE], s = "lambda.min", type = "response"))
  list(prob = prob, selected = sel)
}
run_nested <- function(X, y, patient, reps, shuffle = FALSE) {
  aucs <- numeric(reps); sel_log <- list(); oof <- list()
  for (r in seq_len(reps)) {
    yy <- y
    if (shuffle) for (id in unique(patient)) { i <- which(patient == id); yy[i] <- sample(y[i]) }   # permute WITHIN patient
    f <- patient_folds(patient, n_outer); pred <- rep(NA_real_, length(yy))
    for (j in sort(unique(f))) {
      tr <- f != j; te <- f == j
      fp <- fit_predict(X[tr, , drop = FALSE], yy[tr], patient[tr], X[te, , drop = FALSE])
      pred[te] <- fp$prob; sel_log[[length(sel_log) + 1]] <- fp$selected
    }
    aucs[r] <- as.numeric(auc(roc(yy, pred, levels = c("Normal", "Tumour"), direction = "<", quiet = TRUE)))
    oof[[r]] <- data.frame(rep = r, sample = rownames(X), patient = patient, truth = yy, prob = pred)
    message(sprintf("rep %d/%d  AUC = %.3f", r, reps, aucs[r]))
  }
  list(auc = aucs, sel = sel_log, oof = do.call(rbind, oof))
}

run_layer <- function(name) {
  cfg <- layers[[name]]; cat("\n=========", name, "=========\n")
  set.seed(2026)
  X <- readRDS(cfg$X)
  if (isTRUE(cfg$transpose)) X <- t(X)                                             # samples x features
  meta <- read.csv(cfg$meta, stringsAsFactors = FALSE); meta <- meta[match(rownames(X), meta$sample), ]
  stopifnot(!anyNA(meta$sample))
  y <- factor(meta$class, levels = c("Normal", "Tumour")); patient <- as.character(meta$patient)
  stopifnot(all(table(patient, y) > 0))
  X <- X[, apply(X, 2, function(v) !anyNA(v) && var(v) > 0), drop = FALSE]      # label-free clean-up
  cat("Matrix:", nrow(X), "samples x", ncol(X), "features;", length(unique(patient)), "patients\n")

  real <- run_nested(X, y, patient, n_reps)
  shuf <- run_nested(X, y, patient, n_reps_shuffle, shuffle = TRUE)
  summ <- data.frame(analysis = c("nested patient-grouped CV", "shuffled-label control"),
                     median_AUC = c(median(real$auc), median(shuf$auc)),
                     min_AUC = c(min(real$auc), min(shuf$auc)), max_AUC = c(max(real$auc), max(shuf$auc)),
                     repeats = c(length(real$auc), length(shuf$auc)))
  print(summ)
  write.csv(summ, paste0(cfg$prefix, "_summary.csv"), row.names = FALSE)
  write.csv(real$oof, paste0(cfg$prefix, "_out_of_fold_predictions.csv"), row.names = FALSE)

  n_fits <- length(real$sel); cat("Number of nested selection fits:", n_fits, "\n")
  freq <- sort(table(unlist(real$sel)) / n_fits, decreasing = TRUE)
  fd <- data.frame(feature = names(freq), selection_frequency = as.numeric(freq))
  if (cfg$ensembl) fd$symbol <- as.character(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
        keys = sub("\\.[0-9]+$", "", fd$feature), keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"))
  write.csv(fd, paste0(cfg$prefix, "_selection_frequency.csv"), row.names = FALSE)
  cat("\nTop 20 most frequently selected features:\n"); print(head(fd, 20), row.names = FALSE)

  ## where do the locked panel features stand?
  pn <- read.csv(cfg$panel_file, stringsAsFactors = FALSE)
  panel <- unique(as.character(pn$feature[pn$omics == cfg$omics]))
  hit <- fd$feature %in% panel; if (cfg$ensembl) hit <- hit | fd$symbol %in% panel
  cat("\nLocked-panel", cfg$omics, "features:", length(panel), "\n"); print(fd[hit, ], row.names = FALSE)
  allid <- fd$feature; if (cfg$ensembl) allid <- c(allid, fd$symbol)
  cat("Never selected:", setdiff(panel, allid), "\n")
  invisible(list(summary = summ, freq = fd))
}

results <- lapply(layers_to_run, run_layer); names(results) <- layers_to_run
