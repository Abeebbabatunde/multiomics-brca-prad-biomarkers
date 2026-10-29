# ================================================
# STEP 0.4 — CREATE STARTUP SCRIPT
# This script reloads everything automatically
# every time you open RStudio
# ================================================

startup_code <- '
# ================================================
# 00_STARTUP.R
# RUN THIS EVERY TIME YOU OPEN RSTUDIO
# ================================================

cat("========================================\\n")
cat("   LOADING THESIS PROJECT...\\n")
cat("========================================\\n\\n")

# Set working directory
setwd("MyThesis_MultiOmics")
options(timeout = 7200)  # 2 hour download timeout

# Load all libraries
cat("Loading libraries...\\n")
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(GEOquery)
  library(SummarizedExperiment)
  library(DESeq2)
  library(minfi)
  library(limma)
  library(edgeR)
  library(clusterProfiler)
  library(EnhancedVolcano)
  library(mixOmics)
  library(MOFA2)
  library(ComplexHeatmap)
  library(STRINGdb)
  library(survival)
  library(survminer)
  library(org.Hs.eg.db)
  library(tidyverse)
  library(ggplot2)
  library(pheatmap)
  library(caret)
  library(randomForest)
  library(pROC)
  library(VennDiagram)
})
cat("✅ Libraries loaded\\n\\n")

# Reload all saved files
cat("Reloading saved data...\\n")

load_file <- function(var, path, type="rds", desc=""){
  if(file.exists(path)){
    obj <- if(type=="rds") readRDS(path) else
             read.csv(path, row.names=1,
                      check.names=FALSE)
    assign(var, obj, envir=.GlobalEnv)
    size <- round(file.size(path)/1024/1024, 1)
    cat(sprintf("  ✅ %-35s (%.1f MB)\\n", desc, size))
  } else {
    cat(sprintf("  ⬜ %-35s not yet created\\n", desc))
  }
}

# Breast cancer files
load_file("clinical_brca",
          "01_RawData/Breast/clinical_brca.csv",
          "csv", "Breast clinical data")
load_file("brca_rna",
          "01_RawData/Breast/brca_rna.rds",
          "rds", "Breast RNA-seq object")
load_file("brca_meth",
          "01_RawData/Breast/brca_meth.rds",
          "rds", "Breast methylation object")
load_file("brca_protein",
          "01_RawData/Breast/brca_protein.rds",
          "rds", "Breast proteomics object")

# Prostate cancer files
load_file("clinical_prad",
          "01_RawData/Prostate/clinical_prad.csv",
          "csv", "Prostate clinical data")
load_file("prad_rna",
          "01_RawData/Prostate/prad_rna.rds",
          "rds", "Prostate RNA-seq object")
load_file("prad_meth",
          "01_RawData/Prostate/prad_meth.rds",
          "rds", "Prostate methylation object")

# Processed data
load_file("brca_rna_norm",
          "02_ProcessedData/Breast/brca_rna_normalized.csv",
          "csv", "Breast normalized RNA")
load_file("brca_meth_norm",
          "02_ProcessedData/Breast/brca_meth_normalized.csv",
          "csv", "Breast normalized methylation")
load_file("prad_rna_norm",
          "02_ProcessedData/Prostate/prad_rna_normalized.csv",
          "csv", "Prostate normalized RNA")

# Results
load_file("brca_degs",
          "03_Results/Transcriptomics/brca_DEGs.csv",
          "csv", "Breast DEGs")
load_file("brca_dmps",
          "03_Results/Epigenomics/brca_DMPs.csv",
          "csv", "Breast DMPs")
load_file("prad_degs",
          "03_Results/Transcriptomics/prad_DEGs.csv",
          "csv", "Prostate DEGs")
load_file("prad_dmps",
          "03_Results/Epigenomics/prad_DMPs.csv",
          "csv", "Prostate DMPs")

# Recreate key variables
if(exists("clinical_brca")){
  early_barcodes_brca <- clinical_brca$submitter_id
  cat(sprintf("\\n  ✅ early_barcodes_brca: %d patients\\n",
              length(early_barcodes_brca)))
}
if(exists("clinical_prad")){
  early_barcodes_prad <- clinical_prad$submitter_id
  cat(sprintf("  ✅ early_barcodes_prad: %d patients\\n",
              length(early_barcodes_prad)))
}
if(exists("brca_rna")){
  rna_counts_brca <- assay(brca_rna, "unstranded")
  rna_meta_brca   <- colData(brca_rna)
  cat(sprintf("  ✅ rna_counts_brca: %d genes x %d samples\\n",
              nrow(rna_counts_brca), ncol(rna_counts_brca)))
}
if(exists("brca_meth")){
  meth_matrix_brca <- assay(brca_meth)
  meth_meta_brca   <- colData(brca_meth)
  cat(sprintf("  ✅ meth_matrix_brca: %d probes x %d samples\\n",
              nrow(meth_matrix_brca), ncol(meth_matrix_brca)))
}

# Progress tracker
cat("\\n========================================\\n")
cat("   YOUR RESEARCH PROGRESS\\n")
cat("========================================\\n")

all_steps <- list(
  list("0.2 Packages installed",
       "01_RawData/Breast"),
  list("0.3 Folders created",
       "01_RawData/Breast"),
  list("1.1 Breast clinical data",
       "01_RawData/Breast/clinical_brca.csv"),
  list("1.2 Breast RNA-seq",
       "01_RawData/Breast/brca_rna.rds"),
  list("1.3 Breast methylation",
       "01_RawData/Breast/brca_meth.rds"),
  list("1.4 Prostate clinical data",
       "01_RawData/Prostate/clinical_prad.csv"),
  list("1.5 Prostate RNA-seq",
       "01_RawData/Prostate/prad_rna.rds"),
  list("1.6 Prostate methylation",
       "01_RawData/Prostate/prad_meth.rds"),
  list("2.1 Breast RNA-seq normalized",
       "02_ProcessedData/Breast/brca_rna_normalized.csv"),
  list("2.2 Breast methylation normalized",
       "02_ProcessedData/Breast/brca_meth_normalized.csv"),
  list("2.3 Prostate RNA-seq normalized",
       "02_ProcessedData/Prostate/prad_rna_normalized.csv"),
  list("2.4 Prostate methylation normalized",
       "02_ProcessedData/Prostate/prad_meth_normalized.csv"),
  list("3.1 Breast DEGs",
       "03_Results/Transcriptomics/brca_DEGs.csv"),
  list("3.2 Breast DMPs",
       "03_Results/Epigenomics/brca_DMPs.csv"),
  list("3.3 Prostate DEGs",
       "03_Results/Transcriptomics/prad_DEGs.csv"),
  list("3.4 Prostate DMPs",
       "03_Results/Epigenomics/prad_DMPs.csv"),
  list("4.1 Pathway analysis",
       "03_Results/Transcriptomics/brca_GO_pathways.csv"),
  list("5.1 Venn diagram",
       "04_Figures/VennDiagram_Breast.png"),
  list("5.2 MOFA integration",
       "03_Results/Integration/MOFA_trained.rds"),
  list("5.3 DIABLO model",
       "03_Results/Integration/DIABLO_model.rds"),
  list("6.0 Network analysis",
       "03_Results/Integration/Network_Edges.csv"),
  list("7.0 ML model",
       "03_Results/MachineLearning/RF_model.rds"),
  list("8.0 Validation",
       "03_Results/MachineLearning/Validation_Results.csv"),
  list("9.0 Survival analysis",
       "03_Results/Survival/Survival_Results.csv"),
  list("10.0 Final biomarker panel",
       "05_FinalResults/Final_Biomarker_Panel.csv")
)

done_count <- 0
for(s in all_steps){
  done   <- file.exists(s[[2]])
  status <- ifelse(done, "✅ DONE", "⬜ PENDING")
  cat(sprintf("  %-35s %s\\n", s[[1]], status))
  if(done) done_count <- done_count + 1
}

cat(sprintf("\\n  PROGRESS: %d / %d steps complete\\n",
            done_count, length(all_steps)))
cat("========================================\\n")
cat("  Ready to continue!\\n")
cat("========================================\\n")
'

# Write startup script to file
writeLines(startup_code,
           "MyThesis_MultiOmics/00_Scripts/00_STARTUP.R")

cat("✅ Startup script created!\n")
cat("📁 Location: MyThesis_MultiOmics/00_Scripts/00_STARTUP.R\n")
cat("\nRUN THIS EVERY TIME YOU OPEN RSTUDIO:\n")
cat('source("MyThesis_MultiOmics/00_Scripts/00_STARTUP.R")\n')