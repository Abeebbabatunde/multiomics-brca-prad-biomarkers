# ============================================================
# STEP 4.2 -- PRAD DMP delta-beta recalculation (diagnostic that produced prad_DMPs_corrected_results.*)
# Inputs : 01_RawData/Prostate/prad_meth.rds (RangedSummarizedExperiment, 485,577 x 213: 191 TP + 22 NT)
#          03_Results/Epigenomics/prad_DMPs_full_results.csv  (object `prad_meth`)
# Outputs: 03_Results/Epigenomics/prad_DMPs_corrected_results.csv / .rds
# Console: delta_beta NA 126,865 -> delta_beta_corrected NA 64,426; 62,439 values filled;
#          358,712 overlapping values agree (max |difference| 4.996e-16);
#          the 64,426 probes without a limma P value stay NA.
# NOT SUPPLIED: the line that created `prad_meth$recalc_delta_beta` is shown below as typed in the console
#          (delta_beta_recalc matched by probe_id).
# ============================================================
library(SummarizedExperiment)
prad_raw <- readRDS("01_RawData/Prostate/prad_meth.rds")
meth_mat <- assay(prad_raw)
tumor  <- colData(prad_raw)$shortLetterCode == "TP"
normal <- colData(prad_raw)$shortLetterCode == "NT"
delta_beta_recalc <- rowMeans(meth_mat[, tumor,  drop = FALSE], na.rm = TRUE) -
                     rowMeans(meth_mat[, normal, drop = FALSE], na.rm = TRUE)

prad_meth$recalc_delta_beta <- delta_beta_recalc[match(prad_meth$probe_id, rownames(meth_mat))]

prad_meth$delta_beta_corrected <- prad_meth$delta_beta
idx <- is.na(prad_meth$delta_beta_corrected) & !is.na(prad_meth$recalc_delta_beta)
prad_meth$delta_beta_corrected[idx] <- prad_meth$recalc_delta_beta[idx]
sum(idx)                                        # 62439

write.csv(prad_meth, "03_Results/Epigenomics/prad_DMPs_corrected_results.csv", row.names = FALSE)
saveRDS(prad_meth,   "03_Results/Epigenomics/prad_DMPs_corrected_results.rds")
