# STEP 2.6 — PROSTATE METHYLATION NORMALIZATION (AS RUN). Output: 400,283 probes (as recorded earlier).
# The first attempt failed with "object 'auto_probes' not found"; the definition below was the fix.
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
auto_probes <- rownames(ann)[!ann$chr %in% c("chrX", "chrY")]
cat("Autosomal probes:", length(auto_probes), "\n")

if(!exists("prad_meth")){
  prad_meth        <- readRDS("01_RawData/Prostate/prad_meth.rds")
  meth_matrix_prad <- assay(prad_meth)
  meth_meta_prad   <- colData(prad_meth)
}
miss_p  <- rowMeans(is.na(meth_matrix_prad))
clean_p <- meth_matrix_prad[miss_p < 0.20,]
clean_p <- clean_p[rownames(clean_p) %in% auto_probes,]
prad_meth_norm <- normalizeBetweenArrays(clean_p, method="quantile")
meth_meta_prad$group <- ifelse(meth_meta_prad$sample_type == "Solid Tissue Normal", "Normal","EarlyProstateCancer")
cat("✅ Prostate methylation normalized!\n  Probes:", nrow(prad_meth_norm),"\n")
write.csv(prad_meth_norm, "02_ProcessedData/Prostate/prad_meth_normalized.csv")
saveRDS(as.data.frame(meth_meta_prad), "02_ProcessedData/Prostate/prad_meth_meta.rds")
cat("✅ SAVED: prad_meth_normalized.csv\n")
