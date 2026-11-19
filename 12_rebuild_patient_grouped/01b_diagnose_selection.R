# 01b_diagnose_selection.R
# Diagnostic only (changes nothing). Runs ONE patient subset of the prostate methylation layer and reports,
# stage by stage, how many features each step keeps and how long it takes. Purpose: find out why the real run
# returned 0 features in all 3 subsets.
project_dir <- "C:/Users/USER/Documents/MyThesis_MultiOmics"
setwd(project_dir)
source("rebuild_functions.R")
set.seed(2026)

X <- readRDS("12_rebuild_inputs/X_methylation_PRAD_rebuild.rds")
meta <- read.csv("12_rebuild_inputs/meta_PRAD_methylation_rebuild.csv", stringsAsFactors = FALSE)
meta <- meta[match(rownames(X), meta$sample), ]
y <- factor(meta$class, levels = c("Normal", "Tumour")); patient <- as.character(meta$patient)
X <- X[, apply(X, 2, function(v) !anyNA(v) && var(v) > 0), drop = FALSE]
ids <- unique(patient); pick <- sample(ids, floor(length(ids) * 0.75)); tr <- patient %in% pick
Xtr <- X[tr, , drop = FALSE]; ytr <- y[tr]; ptr <- patient[tr]
cat("Training subset:", nrow(Xtr), "samples,", length(unique(ptr)), "patients;", ncol(Xtr), "features\n")
cat("Class counts:\n"); print(table(ytr))

tm <- function(label, expr) { t0 <- Sys.time(); v <- expr; cat(sprintf("%-40s %6.1f sec\n", label, as.numeric(difftime(Sys.time(), t0, units = "secs")))); v }

keep <- tm("1. paired filter (top 500)", filter_features(Xtr, ytr, ptr, 500))
Xf <- Xtr[, keep, drop = FALSE]
cat("   features kept by filter:", length(keep), "\n")

cv <- tm("2a. cv.glmnet", cv.glmnet(Xf, ytr, family = "binomial", alpha = 1, foldid = patient_folds(ptr, 5)))
cat("   lambda.min =", signif(cv$lambda.min, 3), " non-zero:", sum(as.vector(coef(cv, s = "lambda.min"))[-1] != 0), "\n")
cat("   lambda.1se =", signif(cv$lambda.1se, 3), " non-zero:", sum(as.vector(coef(cv, s = "lambda.1se"))[-1] != 0), "\n")

for (fr in c(0.6, 0.3, 0.1)) {
  s <- tm(sprintf("2b. LASSO stability, freq >= %.1f", fr), sel_lasso(Xf, ytr, ptr, 50, fr))
  cat("   features:", length(s), "\n")
}
for (fr in c(0.6, 0.3, 0.1)) {
  s <- tm(sprintf("3. SVM-RFE stability, freq >= %.1f", fr), sel_svm(Xf, ytr, ptr, B = 20, freq = fr, n_keep = 20))
  cat("   features:", length(s), "\n")
}
sb <- tm("4. Boruta", sel_boruta(Xf, ytr))
cat("   confirmed features:", length(sb), "\n")

cat("\nDistinct warning messages:\n"); print(unique(names(warnings())))
