# 02_breast_diablo_patient_grouped.R
# Step 12.2 (Objective 1): DIABLO for BREAST methylation + RNA (pairs as defined in 00_define_cohorts_and_inputs.R), with folds assigned by PATIENT, a shuffled-label control,
# and the keepX grid EXTENDED DOWNWARD (1, 2, 3, ...) so that an optimum can be found (thesis recommendation 1).
# This mirrors 11_pre_submission_reanalyses/03a_prostate_diablo_patient_grouped.R for the breast cohort.
# Inputs come from 00_define_cohorts_and_inputs.R. Samples are aligned by patient and class across the two layers. Edit project_dir only.

project_dir <- "C:/Users/USER/Documents/MyThesis_MultiOmics"
setwd(project_dir)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("mixOmics", quietly = TRUE)) BiocManager::install("mixOmics", update = FALSE, ask = FALSE)
suppressPackageStartupMessages(library(mixOmics))
set.seed(2026)

n_top_var <- 2000; ncomp <- 2; k_outer <- 4; reps <- 10
keepX_grid <- c(1, 2, 3, 5, 10, 20, 30, 50)     # extended downward; also look at whether the error rises again at the low end

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

## ---- patient-grouped CV with a shuffled-label control ----
patient_folds <- function(patient, k) {
  ids <- unique(patient); k <- min(k, length(ids))
  f <- sample(rep(seq_len(k), length.out = length(ids))); f[match(patient, ids)]
}
diablo_cv <- function(keep_n, shuffle = FALSE) {
  errs <- numeric(reps)
  for (r in seq_len(reps)) {
    yy <- tissue
    if (shuffle) for (id in unique(patient)) { i <- which(patient == id); yy[i] <- sample(tissue[i]) }
    f <- patient_folds(patient, k_outer); wrong <- 0
    for (j in sort(unique(f))) {
      tr <- f != j; te <- f == j
      keepX <- list(methylation = rep(keep_n, ncomp), RNA = rep(keep_n, ncomp))
      fit <- block.splsda(X = lapply(X_list, function(m) m[tr, , drop = FALSE]),
                          Y = droplevels(yy[tr]), ncomp = ncomp, keepX = keepX, design = design)
      p <- predict(fit, newdata = lapply(X_list, function(m) m[te, , drop = FALSE]), dist = "centroids.dist")
      cls <- p$WeightedVote$centroids.dist[, ncomp]
      wrong <- wrong + sum(cls != as.character(yy[te]))
    }
    errs[r] <- wrong / length(yy)
  }
  mean(errs)
}
grid <- do.call(rbind, lapply(keepX_grid, function(k) {
  cat("keepX =", k, "...\n")
  data.frame(keepX = k, error_real = diablo_cv(k), error_shuffled = diablo_cv(k, shuffle = TRUE))
}))
print(grid, row.names = FALSE)
write.csv(grid, "DIABLO_BRCA_patient_grouped_cv.csv", row.names = FALSE)

## ---- final model at the keepX with the lowest real error (smallest keepX if tied) ----
best <- min(grid$keepX[grid$error_real == min(grid$error_real)])
cat("keepX used for the final model:", best, "\n")
final <- block.splsda(X = X_list, Y = tissue, ncomp = ncomp,
                      keepX = list(methylation = rep(best, ncomp), RNA = rep(best, ncomp)), design = design)
sel <- do.call(rbind, lapply(1:ncomp, function(cc) do.call(rbind, lapply(names(X_list), function(b) {
  v <- selectVar(final, block = b, comp = cc)[[b]]$value
  data.frame(component = cc, block = b, feature = rownames(v), loading = v[, 1])
}))))
write.csv(sel, "DIABLO_BRCA_selected_features.csv", row.names = FALSE)
saveRDS(final, "DIABLO_BRCA_final.rds")
