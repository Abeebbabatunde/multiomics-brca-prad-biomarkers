# 04_prostate_methylation_probe_audit.R
# FIX 4: why do some prostate methylation probes have no test result?
# Result: all 64,426 probes without a P value were excluded BEFORE testing (they are absent from the normalised matrix).
# Nothing failed, so no re-run is needed; this script documents the counts for the Methods.
# Inputs: the DMP table (prad_DMPs_full_results.csv; columns probe_id, logFC, delta_beta, P.Value, adj.P.Val, AveExpr, t, B)
#         and the normalised matrix prad_meth_normalized.csv (first column = probe ID).

setwd("C:/Users/USER/Documents/MyThesis_MultiOmics")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in c("minfi", "IlluminaHumanMethylation450kanno.ilmn12.hg19"))
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)
suppressPackageStartupMessages({library(minfi); library(IlluminaHumanMethylation450kanno.ilmn12.hg19)})

full <- read.csv("03_Results/Epigenomics/prad_DMPs_full_results.csv", stringsAsFactors = FALSE)
stopifnot(all(c("probe_id", "adj.P.Val", "delta_beta", "P.Value") %in% names(full)))

f <- "02_ProcessedData/Prostate/prad_meth_normalized.csv"                  # read the probe IDs only (first column)
n <- ncol(read.csv(f, nrows = 1, check.names = FALSE))
ids <- read.csv(f, colClasses = c(NA, rep("NULL", n - 1)), check.names = FALSE)[[1]]

ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
chr <- as.character(ann$chr[match(full$probe_id, rownames(ann))])
snp <- !is.na(ann$Probe_rs[match(full$probe_id, rownames(ann))])
sex <- chr %in% c("chrX", "chrY")

dm <- !is.na(full$adj.P.Val) & full$adj.P.Val < 0.05 & abs(full$delta_beta) >= 0.1
dm[is.na(dm)] <- FALSE
in_mat <- full$probe_id %in% ids
no_P   <- is.na(full$P.Value)

cat("Probes on the array table:", nrow(full), "\n")
cat("Probes in the normalised matrix:", length(ids), "\n")
cat("No P value (not tested):", sum(no_P), " of which absent from the matrix:", sum(no_P & !in_mat), "\n")
cat("Tested:", sum(!no_P), "\n")
cat("Tested but not in the matrix:", sum(!no_P & !in_mat), " (chrX/Y:", sum(!no_P & !in_mat & sex), ")\n")
cat("DMPs (adj.P < 0.05, |delta-beta| >= 0.1):", sum(dm), "\n")
cat("DMPs not in the matrix:", sum(dm & !in_mat), " (chrX/Y:", sum(dm & !in_mat & sex), ")\n")

grp <- ifelse(no_P, "A_no_P_value", ifelse(in_mat, "B_tested_in_matrix", "C_tested_not_in_matrix"))
tab <- data.frame(group = grp, sex_chr = sex, has_SNP = snp, not_cg = !grepl("^cg", full$probe_id))
print(table(grp))
print(aggregate(cbind(sex_chr, has_SNP, not_cg) ~ group, data = tab, FUN = mean))
