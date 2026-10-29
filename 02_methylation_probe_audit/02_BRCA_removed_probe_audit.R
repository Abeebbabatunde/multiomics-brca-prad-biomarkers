# BRCA methylation: what do the 180,024 probes absent from the DMP results correspond to?
# Reconstructed from the successful console run (raw 485,577; results table 305,553; removed 180,024).
# The final NA-rate block is included but its output was NOT seen (first attempt failed: brca_meth_raw is a
# RangedSummarizedExperiment, so assay() is needed). Cross-reactive-probe check was never run.
library(SummarizedExperiment); library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
brca_meth_raw <- readRDS("01_RawData/Breast/brca_meth.rds"); dim(brca_meth_raw)          # 485577 x 649
brca_meth <- readRDS("results/step2_2_methylation_diff/breast_methylation_diffmeth_results.rds")
removed_probes <- setdiff(rownames(brca_meth_raw), brca_meth$Probe); length(removed_probes)  # 180024
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
removed_ann <- ann450k[rownames(ann450k) %in% removed_probes, ]
table(removed_ann$chr %in% c("chrX","chrY"))                       # 4831 sex-chromosome
sum(!is.na(removed_ann$Probe_rs)); sum(!is.na(removed_ann$CpG_rs)); sum(!is.na(removed_ann$SBE_rs))  # 38969 / 15997 / 7391
snp_probes <- rownames(removed_ann)[!is.na(removed_ann$Probe_rs) | !is.na(removed_ann$CpG_rs) | !is.na(removed_ann$SBE_rs)]
length(snp_probes)                                                  # 55141
union_snp_sex <- union(snp_probes, rownames(removed_ann)[removed_ann$chr %in% c("chrX","chrY")])
length(union_snp_sex)                                               # 59267
leftover <- setdiff(removed_probes, union_snp_sex); length(leftover)  # 120757 unexplained
# --- output not seen below ---
lv <- assay(brca_meth_raw[leftover, ])
na_rate <- rowSums(is.na(lv)) / ncol(lv); summary(na_rate)
sum(rowSums(is.na(lv)) > 0); sum(rowSums(is.na(lv)) == ncol(lv))
