# STEP 2.2 — BREAST METHYLATION NORMALIZATION (AS RUN; missing <20%, autosomal only, quantile). Output: 393,747 probes (as recorded earlier).
# NOTE: probes are further filtered (MASK_general / NA) inside 03_differential/03_2_BRCA_methylation_limma_AS_RUN.R.
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

if(!exists("brca_meth")){
  brca_meth        <- readRDS("01_RawData/Breast/brca_meth.rds")
  meth_matrix_brca <- assay(brca_meth)
  meth_meta_brca   <- colData(brca_meth)
}
miss  <- rowMeans(is.na(meth_matrix_brca))
clean <- meth_matrix_brca[miss < 0.20,]
ann         <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
auto_probes <- rownames(ann)[!ann$chr %in% c("chrX","chrY")]
clean       <- clean[rownames(clean) %in% auto_probes,]
brca_meth_norm <- normalizeBetweenArrays(clean, method="quantile")
meth_meta_brca$group <- ifelse(meth_meta_brca$sample_type == "Solid Tissue Normal", "Normal","EarlyBreastCancer")
cat("✅ Probes after QC:", nrow(brca_meth_norm),"\n")
write.csv(brca_meth_norm, "02_ProcessedData/Breast/brca_meth_normalized.csv")
saveRDS(as.data.frame(meth_meta_brca), "02_ProcessedData/Breast/brca_meth_meta.rds")
cat("✅ SAVED: brca_meth_normalized.csv\n")
