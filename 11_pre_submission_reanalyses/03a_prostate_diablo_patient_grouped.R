# 03a_prostate_diablo_patient_grouped.R
# FIX 3 (part 1): DIABLO for prostate methylation + RNA, with folds assigned by PATIENT, the RNA on the log2 scale
# (from 00_build_inputs.R) and a shuffled-label control. Run 03b separately (it needs a fresh R session).
# The 2,000 most variable features per block are chosen without the labels.

setwd("C:/Users/USER/Documents/MyThesis_MultiOmics")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("mixOmics", quietly = TRUE)) BiocManager::install("mixOmics", update = FALSE, ask = FALSE)
suppressPackageStartupMessages(library(mixOmics))
set.seed(2026)

## ---- settings ----
n_top_var <- 2000; ncomp <- 2; k_outer <- 4; reps <- 10
keepX_grid <- c(2, 3, 5, 10, 20, 30, 50)
final_keepX <- 20      # keepX of the reported final model. It had the lowest error in the first grid (2-20), and the error is
                       # flat beyond about 20 (0.044 at 30, 0.024 at 50, within one sample of 34); 20 features per block
                       # keeps the selected list short. Set to NA to use the keepX with the lowest error in the full grid.

## ---- data ----
meta <- read.csv("meta_PRAD_matched.csv", stringsAsFactors = FALSE)
rna  <- readRDS("PRAD_rna_log2.rds")
X    <- readRDS("X_methylation_PRAD.rds")                 # M-values (made in 00_build_inputs.R)
rownames(X) <- substr(rownames(X), 1, 16)
stopifnot(all(meta$sample %in% rownames(X)))
X <- X[meta$sample, ]
X <- X[, apply(X, 2, function(v) !anyNA(v) && var(v) > 0), drop = FALSE]
if (min(X) >= 0 && max(X) <= 1) {                          # safety check: convert if beta values were supplied
  X <- pmin(pmax(X, 0.001), 0.999); X <- log2(X / (1 - X)); cat("Methylation converted from beta to M-values\n")
} else cat("Methylation used as supplied (already M-values)\n")
top_var <- function(m, n) m[, order(apply(m, 2, var), decreasing = TRUE)[seq_len(min(n, ncol(m)))]]
X_list <- list(methylation = top_var(X, n_top_var), RNA = top_var(t(rna[, meta$sample]), n_top_var))
tissue  <- factor(meta$class, levels = c("Normal", "Tumour")); patient <- as.character(meta$patient)
cat("Blocks:", sapply(X_list, ncol), " samples:", nrow(X_list[[1]]), "\n")
design <- matrix(0.1, 2, 2, dimnames = list(names(X_list), names(X_list))); diag(design) <- 0

## ---- patient-grouped cross-validation with a shuffled-label control ----
patient_folds <- function(patient, k) {
  ids <- unique(patient); k <- min(k, length(ids))
  f <- sample(rep(seq_len(k), length.out = length(ids))); f[match(patient, ids)]
}
diablo_cv <- function(keep_n, shuffle = FALSE) {
  errs <- numeric(reps)
  for (r in seq_len(reps)) {
    yy <- tissue
    if (shuffle) for (id in unique(patient)) { i <- which(patient == id); yy[i] <- sample(tissue[i]) }   # within patient
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
write.csv(grid, "DIABLO_PRAD_patient_grouped_cv.csv", row.names = FALSE)

## ---- final model and selected features ----
best <- if (is.na(final_keepX)) grid$keepX[which.min(grid$error_real)] else final_keepX
cat("keepX used for the final model:", best, "\n")
final <- block.splsda(X = X_list, Y = tissue, ncomp = ncomp,
                      keepX = list(methylation = rep(best, ncomp), RNA = rep(best, ncomp)), design = design)
sel <- do.call(rbind, lapply(1:ncomp, function(cc) do.call(rbind, lapply(names(X_list), function(b) {
  v <- selectVar(final, block = b, comp = cc)[[b]]$value
  data.frame(component = cc, block = b, feature = rownames(v), loading = v[, 1])
}))))
write.csv(sel, "DIABLO_PRAD_selected_features.csv", row.names = FALSE)
saveRDS(final, "DIABLO_PRAD_final.rds")

## ---- gene symbols and overlap with the locked prostate panel ----
suppressPackageStartupMessages({library(minfi); library(IlluminaHumanMethylation450kanno.ilmn12.hg19)})
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db", update = FALSE, ask = FALSE)
suppressPackageStartupMessages(library(org.Hs.eg.db))
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
sel$gene <- NA_character_
m <- sel$block == "methylation"; r <- sel$block == "RNA"
g <- ann[sel$feature[m], "UCSC_RefGene_Name"]
sel$gene[m] <- vapply(strsplit(as.character(g), ";"), function(x) paste(unique(x), collapse = ";"), "")
sel$gene[r] <- as.character(mapIds(org.Hs.eg.db, keys = sub("\\..*", "", sel$feature[r]),
                                   column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
sel$gene[sel$gene %in% c("", "NA")] <- NA
panel <- read.csv("10_FinalResults/Step10_FinalPanel/Step10_2_PRAD_final_panel.csv", stringsAsFactors = FALSE)
sel$in_locked_panel <- sel$feature %in% panel$feature | sel$gene %in% panel$feature
print(sel[, c("component", "block", "feature", "gene", "loading", "in_locked_panel")], row.names = FALSE)
cat("Selected features also in the locked panel:", sum(sel$in_locked_panel), "\n")
write.csv(sel, "DIABLO_PRAD_selected_features_annotated.csv", row.names = FALSE)
