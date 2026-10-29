# STEP 3.1 — BREAST DEGs (DESeq2). AS RUN (console d410f2ea): 12,402 DEGs (8,558 up, 3,844 down); outlier replacement refit 8,566 genes.
# STEP 3.1 — BREAST DEGs (Gene Expression)
# Time: 10 to 20 minutes
# ================================================

library(DESeq2)
library(EnhancedVolcano)

if(!exists("dds_brca")){
  dds_brca <- readRDS(
    "02_ProcessedData/Breast/dds_brca.rds")
}

cat("Running DESeq2 for breast cancer...\n")
dds_brca <- DESeq(dds_brca)

res_brca      <- results(dds_brca,
  contrast=c("group","EarlyBreastCancer","Normal"),
  alpha=0.05)
brca_degs_all <- as.data.frame(res_brca)
brca_degs_all <- brca_degs_all[
  complete.cases(brca_degs_all),]

brca_degs <- brca_degs_all[
  brca_degs_all$padj < 0.05 &
  abs(brca_degs_all$log2FoldChange) > 1,]

cat("✅ Breast DEGs:", nrow(brca_degs),"\n")
cat("  Upregulated  :",
    sum(brca_degs$log2FoldChange>0),"\n")
cat("  Downregulated:",
    sum(brca_degs$log2FoldChange<0),"\n")

# Volcano plot
v <- EnhancedVolcano(brca_degs_all,
  lab=rownames(brca_degs_all),
  x="log2FoldChange", y="padj",
  title="Early Breast Cancer — DEGs",
  pCutoff=0.05, FCcutoff=1.0)
print(v)
ggsave("04_Figures/Volcano_Breast_DEGs.pdf",
       v, width=10, height=8)



write.csv(brca_degs,
  "03_Results/Transcriptomics/brca_DEGs.csv")
write.csv(brca_degs_all,
  "03_Results/Transcriptomics/brca_DEGs_all.csv")
cat("✅ SAVED: brca_DEGs.csv\n")
