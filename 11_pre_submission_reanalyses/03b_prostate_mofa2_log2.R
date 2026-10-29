# 03b_prostate_mofa2_log2.R
# FIX 3 (part 2): MOFA2 for prostate methylation + RNA with the RNA on the log2 scale.
# Run in a FRESH R session (Sys.setenv must come before reticulate is loaded).
# MOFA2 needs the Python package mofapy2. On Windows the automatic (basilisk) installation failed here, so Python 3.11
# was installed from python.org and "pip install mofapy2" was run once. Edit the path below to a Python with mofapy2.
# If basilisk works on your machine, delete the Sys.setenv line and set use_basilisk = TRUE.

python_path <- "C:/Users/USER/AppData/Local/Programs/Python/Python311/python.exe"
Sys.setenv(RETICULATE_PYTHON = python_path)
has <- function(p) {
  out <- suppressWarnings(system2(p, c("-c", shQuote("import mofapy2")), stdout = TRUE, stderr = TRUE))
  is.null(attr(out, "status"))
}
if (!has(python_path)) system2(python_path, c("-m", "pip", "install", "mofapy2"))
stopifnot("mofapy2 is not available in this Python" = has(python_path))

setwd("C:/Users/USER/Documents/MyThesis_MultiOmics")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("MOFA2", quietly = TRUE)) BiocManager::install("MOFA2", update = FALSE, ask = FALSE)
suppressPackageStartupMessages(library(MOFA2))
set.seed(2026)

## ---- same data preparation as DIABLO ----
meta <- read.csv("meta_PRAD_matched.csv", stringsAsFactors = FALSE)
rna  <- readRDS("PRAD_rna_log2.rds")
X    <- readRDS("X_methylation_PRAD.rds")
rownames(X) <- substr(rownames(X), 1, 16)
X <- X[meta$sample, ]
X <- X[, apply(X, 2, function(v) !anyNA(v) && var(v) > 0), drop = FALSE]
if (min(X) >= 0 && max(X) <= 1) { X <- pmin(pmax(X, 0.001), 0.999); X <- log2(X / (1 - X)) }
top_var <- function(m, n) m[, order(apply(m, 2, var), decreasing = TRUE)[seq_len(min(n, ncol(m)))]]
meth_top <- top_var(X, 2000); rna_top <- top_var(t(rna[, meta$sample]), 2000)
tissue  <- factor(meta$class, levels = c("Normal", "Tumour")); patient <- as.character(meta$patient)

## ---- MOFA2 ----
mofa <- create_mofa(list(methylation = t(meth_top), RNA = t(rna_top)))
data_opts  <- get_default_data_options(mofa);  data_opts$scale_views <- TRUE
model_opts <- get_default_model_options(mofa); model_opts$num_factors <- 5
train_opts <- get_default_training_options(mofa); train_opts$seed <- 2026
mofa <- prepare_mofa(mofa, data_options = data_opts, model_options = model_opts, training_options = train_opts)
mofa <- run_mofa(mofa, outfile = "MOFA_PRAD_log2.hdf5", use_basilisk = FALSE)

## ---- variance explained and tumour-versus-normal test for each factor ----
r2 <- calculate_variance_explained(mofa)
cat("\nTotal variance explained per view (%):\n"); print(round(r2$r2_total[[1]], 1))
cat("\nVariance explained per factor and view (%):\n"); print(round(r2$r2_per_factor[[1]], 1))
Z <- get_factors(mofa)[[1]]; Z <- Z[meta$sample, , drop = FALSE]
mean_rna <- rowMeans(rna_top)
paired_p <- function(z) {
  zt <- tapply(z[tissue == "Tumour"], patient[tissue == "Tumour"], mean)
  zn <- tapply(z[tissue == "Normal"], patient[tissue == "Normal"], mean)
  pp <- intersect(names(zt), names(zn))
  wilcox.test(zt[pp], zn[pp], paired = TRUE, exact = FALSE)$p.value
}
cat("\nPaired tumour vs normal, per factor (Wilcoxon signed-rank):\n")
for (k in seq_len(ncol(Z)))
  cat(sprintf("Factor %d: P = %.2g   (Spearman with mean RNA signal: %.2f)\n", k, paired_p(Z[, k]),
              cor(Z[, k], mean_rna, method = "spearman")))
cat("Factor 1 after adjusting for the mean RNA signal: P =", signif(paired_p(resid(lm(Z[, 1] ~ mean_rna))), 2), "\n")
saveRDS(Z, "MOFA_PRAD_factors.rds")

## ---- strongest RNA weights of Factors 1 and 2, with gene symbols ----
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db", update = FALSE, ask = FALSE)
suppressPackageStartupMessages(library(org.Hs.eg.db))
W <- get_weights(mofa, views = "RNA", factors = 1:2)[[1]]
sym <- function(ids) as.character(mapIds(org.Hs.eg.db, keys = sub("\\..*", "", ids),
                                         column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
for (k in 1:2) {
  o <- order(abs(W[, k]), decreasing = TRUE)[1:15]
  cat("\nFactor", k, "- top 15 RNA genes by absolute weight\n")
  print(data.frame(gene = sym(rownames(W)[o]), id = rownames(W)[o], weight = round(W[o, k], 2)), row.names = FALSE)
}
