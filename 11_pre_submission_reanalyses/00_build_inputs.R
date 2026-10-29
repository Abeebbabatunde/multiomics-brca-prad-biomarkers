# 00_build_inputs.R
# Builds the input files used by scripts 01 to 04:
#   X_methylation_PRAD.rds + meta_PRAD.csv            (prostate methylation, M-values, 34 samples)
#   X_methylation_BRCA.rds + meta_BRCA.csv            (breast methylation, M-values, 102 samples)
#   X_rna_BRCA.rds + meta_BRCA_rna.csv                (breast RNA, variance-stabilised, 102 samples)
#   PRAD_rna_log2.rds + meta_PRAD_matched.csv         (prostate RNA, log2(normalised count + 1), 34 samples)
# Run once, from the project folder. Edit project_dir only.

project_dir <- "C:/Users/USER/Documents/MyThesis_MultiOmics"
setwd(project_dir)
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
library(data.table)

to_M        <- function(beta) { beta <- pmin(pmax(beta, 0.001), 0.999); log2(beta / (1 - beta)) }
label_class <- function(g) ifelse(grepl("Normal", g, ignore.case = TRUE), "Normal", "Tumour")

## ---------------- 1. Prostate methylation ----------------
man <- read.csv("05_FinalResults/Step5.0_Integration/PRAD/cohort_meth_rna.csv")
man <- man[!is.na(man$meth), ]
m <- fread("02_ProcessedData/Prostate/prad_meth_normalized.csv", data.table = FALSE)
rownames(m) <- m[[1]]; m <- as.matrix(m[, -1])                 # probes x samples (beta values)
meta <- data.frame(sample = man$meth, patient = man$patient, class = label_class(man$group))
cat("Prostate: manifest samples found in matrix:", sum(meta$sample %in% colnames(m)), "of", nrow(meta), "\n")
stopifnot(all(meta$sample %in% colnames(m)))
m34 <- m[, meta$sample]
m34 <- m34[rowSums(is.na(m34)) == 0, ]                         # probes with complete data in the 34 paired samples
X <- t(to_M(m34))                                              # samples x probes, M-values
saveRDS(X, "X_methylation_PRAD.rds"); write.csv(meta, "meta_PRAD.csv", row.names = FALSE)
print(table(meta$class)); print(dim(X))

## ---------------- 2. Breast methylation ----------------
man_b <- read.csv("05_FinalResults/Step5.0_Integration/BRCA/cohort_meth_rna.csv")
man_b <- man_b[!is.na(man_b$meth), ]
mb <- fread("02_ProcessedData/Breast/brca_meth_normalized.csv", data.table = FALSE)
barcode_cols <- grepl("^TCGA", names(mb)); id_col <- which(!barcode_cols)[1]
rownames(mb) <- mb[[id_col]]; mb <- as.matrix(mb[, barcode_cols])
meta_b <- data.frame(sample = man_b$meth, patient = man_b$patient, class = label_class(man_b$group))
cat("Breast methylation: manifest samples found in matrix:", sum(meta_b$sample %in% colnames(mb)), "of", nrow(meta_b), "\n")
stopifnot(all(meta_b$sample %in% colnames(mb)))
mb34 <- mb[, meta_b$sample]
mb34 <- mb34[rowSums(is.na(mb34)) == 0, ]
Xb <- t(to_M(mb34))
saveRDS(Xb, "X_methylation_BRCA.rds"); write.csv(meta_b, "meta_BRCA.csv", row.names = FALSE)
print(table(meta_b$class)); print(dim(Xb))

## ---------------- 3. Breast RNA (already variance-stabilised) ----------------
v <- readRDS("02_ProcessedData/Breast/brca_vst.rds")
if (inherits(v, "SummarizedExperiment")) v <- SummarizedExperiment::assay(v)
v <- as.matrix(v)                                              # genes x samples
man_r <- read.csv("05_FinalResults/Step5.0_Integration/BRCA/cohort_meth_rna.csv")
man_r <- man_r[!is.na(man_r$rna), ]
meta_r <- data.frame(sample = man_r$rna, patient = man_r$patient, class = label_class(man_r$group))
cat("Breast RNA: manifest samples found in matrix:", sum(meta_r$sample %in% colnames(v)), "of", nrow(meta_r), "\n")
stopifnot(all(meta_r$sample %in% colnames(v)))
vv <- v[, meta_r$sample]
vv <- vv[rowSums(is.na(vv)) == 0, ]
Xr <- t(vv)
saveRDS(Xr, "X_rna_BRCA.rds"); write.csv(meta_r, "meta_BRCA_rna.csv", row.names = FALSE)
print(table(meta_r$class)); print(dim(Xr))

## ---------------- 4. Prostate RNA on a log2 scale (fix 3) ----------------
# The file holds normalised counts (not raw counts), so log2(count + 1) is used instead of vst.
r <- read.csv("02_ProcessedData/Prostate/prad_rna_normalized.csv", check.names = FALSE, stringsAsFactors = FALSE)
genes <- as.character(r[[1]])
rm0 <- as.matrix(r[, -1]); mode(rm0) <- "numeric"; rownames(rm0) <- make.unique(genes)
meta_p <- read.csv("meta_PRAD.csv", stringsAsFactors = FALSE)
if (sum(colnames(rm0) %in% meta_p$sample) < nrow(meta_p)) {     # barcodes match on the first 16 characters
  colnames(rm0) <- substr(colnames(rm0), 1, 16)
  meta_p$sample <- substr(meta_p$sample, 1, 16)
}
stopifnot(all(meta_p$sample %in% colnames(rm0)))
rm0 <- rm0[, meta_p$sample]
keep <- rowSums(rm0 >= 10) >= 5                                 # same low-count rule as breast: >= 10 in >= 5 samples
rna_log <- log2(rm0[keep, ] + 1)
sp <- function(m) round(cor(rowMeans(m), apply(m, 1, sd), method = "spearman"), 2)
cat("Genes kept:", sum(keep), "of", nrow(rm0), "\n")
cat("Mean-SD Spearman correlation: original", sp(rm0[keep, ]), " log2", sp(rna_log), "\n")
saveRDS(rna_log, "PRAD_rna_log2.rds"); write.csv(meta_p, "meta_PRAD_matched.csv", row.names = FALSE)
cat("Inputs written to", getwd(), "\n")
