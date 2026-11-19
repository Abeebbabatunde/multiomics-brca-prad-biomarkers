# 02b_leakage_demonstration.R
# Step 12.2b (Objective 1): the SAME breast DIABLO model on the SAME rebuilt data, cross-validated two ways:
#   (a) patient-grouped folds : a patient's tumour and normal samples always stay together (correct);
#   (b) sample-level folds    : samples are assigned to folds at random, so a test sample's matched partner (same patient) is often in the training set (leaky).
# Everything else is identical (data, 2000 most variable features per layer, ncomp, design, 4 folds x 10 repeats, keepX), so the difference in error
# is the effect of patient leakage alone. A shuffled-label control (labels permuted within patient) is run under both fold types.
# Inputs come from 00_define_cohorts_and_inputs.R. Edit project_dir only.


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


## ---- two ways of making folds ----
patient_folds <- function(patient, k) {
  ids <- unique(patient); k <- min(k, length(ids))
  f <- sample(rep(seq_len(k), length.out = length(ids))); f[match(patient, ids)]
}
sample_folds <- function(n, k) sample(rep(seq_len(k), length.out = n))              # ignores patients
partner_apart <- function(f) mean(sapply(seq_along(f), function(i) any(f[patient == patient[i] & seq_along(f) != i] != f[i])))
# = share of samples whose matched partner sits in a different fold (so is in the training set when the sample is tested)

cv_run <- function(keep_n, leaky, shuffle = FALSE) {
  errs <- numeric(reps); share <- numeric(reps)
  for (r in seq_len(reps)) {
    yy <- tissue
    if (shuffle) for (id in unique(patient)) { i <- which(patient == id); yy[i] <- sample(tissue[i]) }
    f <- if (leaky) sample_folds(length(yy), k_outer) else patient_folds(patient, k_outer)
    share[r] <- partner_apart(f); wrong <- 0
    for (j in sort(unique(f))) {
      tr <- f != j; te <- f == j
      keepX <- list(methylation = rep(keep_n, ncomp), RNA = rep(keep_n, ncomp))
      fit <- block.splsda(X = lapply(X_list, function(m) m[tr, , drop = FALSE]),
                          Y = droplevels(yy[tr]), ncomp = ncomp, keepX = keepX, design = design)
      p <- predict(fit, newdata = lapply(X_list, function(m) m[te, , drop = FALSE]), dist = "centroids.dist")
      wrong <- wrong + sum(p$WeightedVote$centroids.dist[, ncomp] != as.character(yy[te]))
    }
    errs[r] <- wrong / length(yy)
  }
  c(mean = mean(errs), sd = sd(errs), partner_in_training = mean(share))
}

keepX_demo <- c(1, 5, 20, 50)
res <- do.call(rbind, lapply(keepX_demo, function(k) {
  cat("keepX =", k, "...\n")
  a <- cv_run(k, leaky = FALSE); b <- cv_run(k, leaky = TRUE)
  a0 <- cv_run(k, leaky = FALSE, shuffle = TRUE); b0 <- cv_run(k, leaky = TRUE, shuffle = TRUE)
  data.frame(keepX = k,
             error_patient_grouped = a["mean"], sd_patient_grouped = a["sd"],
             error_sample_level = b["mean"], sd_sample_level = b["sd"],
             share_partner_in_training_sample_level = b["partner_in_training"],
             shuffled_patient_grouped = a0["mean"], shuffled_sample_level = b0["mean"], row.names = NULL)
}))
print(res, row.names = FALSE)
write.csv(res, "DIABLO_BRCA_leakage_demonstration.csv", row.names = FALSE)
cat("\nLeakage effect = error_patient_grouped minus error_sample_level (same data, same model). share_partner_in_training_sample_level is the proportion of\n",
    "test samples whose matched tumour/normal partner was in the training set under sample-level folds (0 under patient-grouped folds by construction).\n", sep = "")
