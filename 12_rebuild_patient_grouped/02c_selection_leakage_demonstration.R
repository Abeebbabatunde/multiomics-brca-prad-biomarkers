# 02c_selection_leakage_demonstration.R
# Step 12.2c (Objective 1): where could an over-optimistic error come from? Here ALL folds are patient-grouped (correct), and only the
# timing of the LABEL-USING feature filter differs:
#   INSIDE  : the paired moderated-t filter (top n_filter features per layer, patient as blocking factor) is recomputed on the TRAINING patients of each fold only (correct);
#   OUTSIDE : the same filter is computed ONCE on ALL samples, labels included, before cross-validation (selection leakage).
# Then the same DIABLO model is cross-validated on the selected features. A shuffled-label control (labels permuted within patient) is run for both:
# a sound procedure stays near 50% error on shuffled labels; selection leakage pushes the error below 50% even though the labels carry no signal.
# Inputs come from 00_define_cohorts_and_inputs.R. Edit project_dir only.

project_dir <- "C:/Users/USER/Documents/MyThesis_MultiOmics"
setwd(project_dir)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("mixOmics", quietly = TRUE)) BiocManager::install("mixOmics", update = FALSE, ask = FALSE)
suppressPackageStartupMessages({library(mixOmics); library(limma)})
set.seed(2026)

n_top_var <- 2000; ncomp <- 2; k_outer <- 4; reps <- 10
n_prefilter <- 20000     # label-free pre-reduction (most variable features per layer) so that the supervised filter runs quickly
n_filter    <- 500       # features kept per layer by the supervised (paired moderated-t) filter

## ---- data: align methylation and RNA by patient and class ----
Xm <- readRDS("12_rebuild_inputs/X_methylation_BRCA_rebuild.rds"); mm <- read.csv("12_rebuild_inputs/meta_BRCA_methylation_rebuild.csv", stringsAsFactors = FALSE)
Xr <- readRDS("12_rebuild_inputs/X_rna_BRCA_rebuild.rds");       mr <- read.csv("12_rebuild_inputs/meta_BRCA_rna_rebuild.csv", stringsAsFactors = FALSE)
mm <- mm[match(rownames(Xm), mm$sample), ]; mr <- mr[match(rownames(Xr), mr$sample), ]
mm$key <- paste(mm$patient, mm$class); mr$key <- paste(mr$patient, mr$class)
mm <- mm[!duplicated(mm$key), ]; mr <- mr[!duplicated(mr$key), ]            # one sample per patient and class
keys <- intersect(mm$key, mr$key)
mm <- mm[match(keys, mm$key), ]; mr <- mr[match(keys, mr$key), ]
Xm <- Xm[mm$sample, ]; Xr <- Xr[mr$sample, ]
Xm <- Xm[, apply(Xm, 2, function(v) !anyNA(v) && var(v) > 0), drop = FALSE]
Xr <- Xr[, apply(Xr, 2, function(v) !anyNA(v) && var(v) > 0), drop = FALSE]
top_var <- function(m, n) m[, order(apply(m, 2, var), decreasing = TRUE)[seq_len(min(n, ncol(m)))]]   # label-free
X_list <- list(methylation = top_var(Xm, n_top_var), RNA = top_var(Xr, n_top_var))
rownames(X_list$methylation) <- rownames(X_list$RNA) <- keys
tissue <- factor(mm$class, levels = c("Normal", "Tumour")); patient <- as.character(mm$patient)
cat("Samples:", nrow(X_list[[1]]), "; patients:", length(unique(patient)), "; blocks:", sapply(X_list, ncol), "\n")
stopifnot(all(table(patient, tissue) > 0))
design <- matrix(0.1, 2, 2, dimnames = list(names(X_list), names(X_list))); diag(design) <- 0



## ---- label-free pre-reduction, then the supervised filter ----
X_pool <- list(methylation = top_var(Xm, n_prefilter), RNA = top_var(Xr, n_prefilter))
rownames(X_pool$methylation) <- rownames(X_pool$RNA) <- keys
filter_features <- function(Xtr, ytr, ptr, n_keep) {          # paired moderated t, patient as blocking factor (same as rebuild_functions.R)
  design_f <- model.matrix(~ factor(ptr) + ytr)
  fit <- eBayes(lmFit(t(Xtr), design_f))
  tt <- topTable(fit, coef = ncol(design_f), number = Inf, sort.by = "P")
  rownames(tt)[seq_len(min(n_keep, nrow(tt)))]
}
get_feats <- function(rows, yy) lapply(X_pool, function(m) filter_features(m[rows, , drop = FALSE], droplevels(yy[rows]), patient[rows], n_filter))
patient_folds <- function(patient, k) {
  ids <- unique(patient); k <- min(k, length(ids))
  f <- sample(rep(seq_len(k), length.out = length(ids))); f[match(patient, ids)]
}
all_rows <- rep(TRUE, length(tissue))

cv_sel <- function(keep_n, outside, shuffle = FALSE) {
  errs <- numeric(reps)
  for (r in seq_len(reps)) {
    yy <- tissue
    if (shuffle) for (id in unique(patient)) { i <- which(patient == id); yy[i] <- sample(tissue[i]) }
    f <- patient_folds(patient, k_outer); wrong <- 0
    feats_all <- if (outside) get_feats(all_rows, yy) else NULL          # labels of ALL samples used here when outside = TRUE
    for (j in sort(unique(f))) {
      tr <- f != j; te <- f == j
      ft <- if (outside) feats_all else get_feats(tr, yy)                # inside: training patients only
      Xtr <- lapply(names(X_pool), function(b) X_pool[[b]][tr, ft[[b]], drop = FALSE]); names(Xtr) <- names(X_pool)
      Xte <- lapply(names(X_pool), function(b) X_pool[[b]][te, ft[[b]], drop = FALSE]); names(Xte) <- names(X_pool)
      keepX <- list(methylation = rep(keep_n, ncomp), RNA = rep(keep_n, ncomp))
      fit <- block.splsda(X = Xtr, Y = droplevels(yy[tr]), ncomp = ncomp, keepX = keepX, design = design)
      p <- predict(fit, newdata = Xte, dist = "centroids.dist")
      wrong <- wrong + sum(p$WeightedVote$centroids.dist[, ncomp] != as.character(yy[te]))
    }
    errs[r] <- wrong / length(yy)
  }
  c(mean = mean(errs), sd = sd(errs))
}

keepX_demo <- c(5, 20)
res <- do.call(rbind, lapply(keepX_demo, function(k) {
  cat("keepX =", k, "...\n")
  a <- cv_sel(k, outside = FALSE);  b <- cv_sel(k, outside = TRUE)
  a0 <- cv_sel(k, outside = FALSE, shuffle = TRUE); b0 <- cv_sel(k, outside = TRUE, shuffle = TRUE)
  data.frame(keepX = k,
             error_filter_inside = a["mean"], sd_inside = a["sd"],
             error_filter_outside = b["mean"], sd_outside = b["sd"],
             shuffled_filter_inside = a0["mean"], shuffled_filter_outside = b0["mean"], row.names = NULL)
}))
print(res, row.names = FALSE)
write.csv(res, "DIABLO_BRCA_selection_leakage_demonstration.csv", row.names = FALSE)
cat("\nInside = filter learned on training patients only (correct). Outside = filter learned on all patients before cross-validation (leaky).\n",
    "On shuffled labels the correct procedure should stay near 0.50; a clearly lower value for 'outside' is selection leakage.\n", sep = "")
