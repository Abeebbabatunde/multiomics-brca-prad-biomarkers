# STEP 2.5 — PROSTATE RNA NORMALIZATION (AS RUN). Size-factor-normalised counts (count scale), genes with counts>0 in >=10% of samples.
library(DESeq2)
if(!exists("prad_rna")){
  prad_rna        <- readRDS("01_RawData/Prostate/prad_rna.rds")
  rna_counts_prad <- assay(prad_rna,"unstranded")
  rna_meta_prad   <- colData(prad_rna)
}
rna_meta_prad$group <- ifelse(rna_meta_prad$sample_type == "Solid Tissue Normal", "Normal","EarlyProstateCancer")
keep_p    <- rowSums(rna_counts_prad>0) >= round(ncol(rna_counts_prad)*0.10)
prad_filt <- rna_counts_prad[keep_p,]
sample_p  <- data.frame(group=factor(rna_meta_prad$group, levels=c("Normal","EarlyProstateCancer")), row.names=colnames(prad_filt))
dds_prad      <- DESeqDataSetFromMatrix(countData=prad_filt, colData=sample_p, design=~group)
dds_prad      <- estimateSizeFactors(dds_prad)
prad_rna_norm <- counts(dds_prad,normalized=TRUE)
cat("✅ Prostate RNA-seq normalized!\n")
saveRDS(dds_prad, "02_ProcessedData/Prostate/dds_prad.rds")
write.csv(prad_rna_norm, "02_ProcessedData/Prostate/prad_rna_normalized.csv")
saveRDS(as.data.frame(rna_meta_prad), "02_ProcessedData/Prostate/prad_rna_meta.rds")
cat("✅ SAVED: prad_rna_normalized.csv\n")
# CAUTION: every non-"Solid Tissue Normal" sample is labelled EarlyProstateCancer; group counts of the 219 samples were not shown.
