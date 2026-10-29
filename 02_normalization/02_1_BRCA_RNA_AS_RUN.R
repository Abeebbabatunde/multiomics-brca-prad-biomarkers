# STEP 2.1 — BREAST RNA NORMALIZATION (AS RUN; assembled from the pasted block + the pasted save lines)
library(DESeq2)
library(ggplot2)

if(!exists("brca_rna")){
  brca_rna        <- readRDS(
    "01_RawData/Breast/brca_rna.rds")
  rna_counts_brca <- assay(brca_rna,"unstranded")
  rna_meta_brca   <- colData(brca_rna)
}

# Group labels
rna_meta_brca$group <- ifelse(
  rna_meta_brca$sample_type ==
  "Solid Tissue Normal",
  "Normal","EarlyBreastCancer"
)

cat("Cancer samples:",
    sum(rna_meta_brca$group=="EarlyBreastCancer"),"\n")
cat("Normal samples:",
    sum(rna_meta_brca$group=="Normal"),"\n")

# Filter low expressed genes
keep <- rowSums(rna_counts_brca > 0) >=
        round(ncol(rna_counts_brca)*0.10)
rna_filt <- rna_counts_brca[keep,]

# DESeq2 normalization
sample_df <- data.frame(
  group = factor(rna_meta_brca$group,
    levels=c("Normal","EarlyBreastCancer")),
  row.names=colnames(rna_filt)
)

dds_brca      <- DESeqDataSetFromMatrix(
  countData=rna_filt,
  colData=sample_df,
  design=~group)
dds_brca      <- estimateSizeFactors(dds_brca)
brca_rna_norm <- counts(dds_brca, normalized=TRUE)

# PCA check
vsd <- vst(dds_brca, blind=FALSE)
p   <- plotPCA(vsd, intgroup="group",
               returnData=TRUE)
pct <- round(100*attr(p,"percentVar"))

pca_plot <- ggplot(p, aes(PC1,PC2,color=group))+
  geom_point(size=3,alpha=0.7)+
  xlab(paste0("PC1: ",pct[1],"% var"))+
  ylab(paste0("PC2: ",pct[2],"% var"))+
  scale_color_manual(values=c(
    "Normal"="#2ECC71",
    "EarlyBreastCancer"="#E74C3C"))+
  theme_bw()+
  ggtitle("PCA: Breast RNA-seq")
print(pca_plot)
ggsave("04_Figures/PCA_Breast_RNA.pdf",
       pca_plot,width=8,height=6)

cat("✅ Breast RNA-seq normalized!\n")

saveRDS(vsd, "02_ProcessedData/Breast/brca_vst.rds")
saveRDS(pca_plot, "04_Figures/PCA_Breast_RNA.rds")
saveRDS(dds_brca, "02_ProcessedData/Breast/dds_brca.rds")
write.csv(brca_rna_norm, "02_ProcessedData/Breast/brca_rna_normalized.csv")
saveRDS(as.data.frame(rna_meta_brca), "02_ProcessedData/Breast/brca_rna_meta.rds")
cat("✅ SAVED: brca_rna_normalized.csv\n")
# CAUTION: group = "Normal" only for sample_type == "Solid Tissue Normal"; every other sample type (e.g. Metastatic) is labelled EarlyBreastCancer.
# The sample_type composition of the 905 BRCA RNA columns was not shown.
