# Step 0.2 package check and Step 0.3 folder check, as pasted (console: Bioconductor 3.23, BiocManager 1.30.27, R 4.6.0 (2026-04-24 ucrt)).
all_pkgs <- c("GEOquery","TCGAbiolinks","minfi","DESeq2","limma","edgeR","clusterProfiler","EnhancedVolcano",
  "mixOmics","MOFA2","ComplexHeatmap","STRINGdb","survminer","survival","org.Hs.eg.db","AnnotationDbi",
  "SummarizedExperiment","tidyverse","ggplot2","pheatmap","caret","randomForest","pROC","VennDiagram","e1071","glmnet")
missing <- c()
for(pkg in all_pkgs){ ok <- requireNamespace(pkg, quietly=TRUE); if(!ok) missing <- c(missing, pkg) }
# If missing: BiocManager::install(missing, ask=FALSE)
# ComplexHeatmap / STRINGdb were re-checked with: sapply(c("ComplexHeatmap","STRINGdb"), requireNamespace, quietly=TRUE)

required <- c("01_RawData/Breast","01_RawData/Prostate","02_ProcessedData/Breast","02_ProcessedData/Prostate",
  "03_Results/Epigenomics","03_Results/Transcriptomics","03_Results/Proteomics","03_Results/Metabolomics",
  "03_Results/Integration","03_Results/MachineLearning","03_Results/Survival","04_Figures","05_FinalResults")
# (the dir.create step 0.3 that created these was not pasted; this is the check)
for(f in required) cat(f, dir.exists(f), "\n")
