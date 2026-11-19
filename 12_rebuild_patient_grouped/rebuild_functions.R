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
