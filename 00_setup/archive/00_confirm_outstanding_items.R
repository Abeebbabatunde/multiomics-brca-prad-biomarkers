# ============================================================
# CONFIRM OUTSTANDING ITEMS  (read-only: changes no data files)
# Run from the project root:  setwd("C:/Users/USER/Documents/MyThesis_MultiOmics")
# Output goes to the console AND to 00_Scripts/confirm_outstanding_output.txt  -> send me that file.
# ============================================================
suppressPackageStartupMessages({library(SummarizedExperiment)})
dir.create("00_Scripts", showWarnings = FALSE)
sink("00_Scripts/confirm_outstanding_output.txt", split = TRUE)
safe <- function(label, expr){
  cat("\n------------------------------------------------\n", label, "\n", sep = "")
  tryCatch(expr, error = function(e) cat("  !! could not run:", conditionMessage(e), "\n"))
}
info <- function(f) if(file.exists(f)) cat(sprintf("  %-60s %8.1f MB  modified %s\n", f,
          file.size(f)/1024/1024, format(file.mtime(f), "%Y-%m-%d %H:%M"))) else cat("  MISSING:", f, "\n")
cat("R:", R.version.string, "\nDate:", format(Sys.time()), "\nWD:", getwd(), "\n")

# 1. SAMPLE-TYPE COMPOSITION (RNA / methylation) ---------------------------------
safe("1a. BRCA RNA (brca_rna.rds): sample_type", {
  x <- readRDS("01_RawData/Breast/brca_rna.rds"); cat("  dim:", dim(x), "\n")
  print(table(colData(x)$sample_type, useNA = "ifany"))
  cat("  patients:", length(unique(substr(colnames(x), 1, 12))), "\n")
})
safe("1b. BRCA methylation (brca_meth.rds): sample_type", {
  x <- readRDS("01_RawData/Breast/brca_meth.rds"); cat("  dim:", dim(x), "\n")
  print(table(colData(x)$sample_type, useNA = "ifany"))
})
safe("1c. PRAD RNA (prad_rna.rds): sample_type", {
  x <- readRDS("01_RawData/Prostate/prad_rna.rds"); cat("  dim:", dim(x), "\n")
  print(table(colData(x)$sample_type, useNA = "ifany"))
})
safe("1d. PRAD methylation (prad_meth.rds): sample_type", {
  x <- readRDS("01_RawData/Prostate/prad_meth.rds"); cat("  dim:", dim(x), "\n")
  print(table(colData(x)$sample_type, useNA = "ifany"))
})
safe("1e. Group tables actually used downstream (02_ProcessedData meta files)", {
  for(f in c("02_ProcessedData/Breast/brca_rna_meta.rds","02_ProcessedData/Breast/brca_meth_meta.rds",
             "02_ProcessedData/Prostate/prad_rna_meta.rds","02_ProcessedData/Prostate/prad_meth_meta.rds")){
    cat(" ", f, "\n"); if(file.exists(f)){ m <- readRDS(f); print(table(m$sample_type, m$group, useNA="ifany")) } else cat("   missing\n")
  }
})

# 2. BRCA RPPA FILE brca_protein_raw.csv (what is it, and how was it built?) -----
safe("2. brca_protein_raw.csv", {
  f <- "01_RawData/Breast/brca_protein_raw.csv"; info(f)
  cat("  first 3 lines:\n"); print(readLines(f, n = 3))
  d <- read.csv(f, stringsAsFactors = FALSE); cat("  dim:", dim(d), "\n  columns:", paste(colnames(d), collapse=", "), "\n")
  if(all(c("sample_barcode","sample_type","peptide_target") %in% colnames(d))){
    cat("  antibodies:", length(unique(d$peptide_target)), " samples:", length(unique(d$sample_barcode)), "\n")
    print(table(d$sample_type[!duplicated(d$sample_barcode)]))
    pt <- unique(d[, c("sample_barcode","sample_type")]); pt$patient <- substr(pt$sample_barcode, 1, 12)
    tab <- table(pt$patient, pt$sample_type)
    cat("  patients with both Primary Tumor and Solid Tissue Normal:", sum(rowSums(tab > 0) == 2), "\n")
    cat("  NA values:", sum(is.na(d$value)), "\n")
  }
})
safe("2b. Proteome Profiling stage of the 30 RPPA patients (needs clinical_brca.csv)", {
  d <- read.csv("01_RawData/Breast/brca_protein_raw.csv"); cl <- read.csv("01_RawData/Breast/clinical_brca.csv")
  p <- unique(substr(d$sample_barcode, 1, 12)); print(table(cl$ajcc_pathologic_stage[cl$submitter_id %in% p], useNA="ifany"))
  cat("  patients not in clinical_brca.csv:", sum(!p %in% cl$submitter_id), "\n")
})

# 3. iTRAQ FILTER RESULT -----------------------------------------------------------
safe("3. iTRAQ brca_protein_early.csv", {
  f <- "01_RawData/Breast/brca_protein_early.csv"; info(f)
  d <- read.csv(f, row.names = 1, check.names = FALSE); cat("  proteins x samples:", dim(d), "\n")
  s <- gsub("\\.", "-", colnames(d)); cat("  sample-type codes (chars 14-15):\n"); print(table(substr(s, 14, 15), useNA="ifany"))
  cl <- read.csv("01_RawData/Breast/clinical_brca.csv")
  m <- cl[cl$submitter_id %in% substr(s, 1, 12), ]; cat("  matched patients:", nrow(m), "\n"); print(table(m$ajcc_pathologic_stage, useNA="ifany"))
  cat("  first sample names:", head(colnames(d), 4), "\n")
})
safe("3b. Raw .cct file present?", {
  print(list.files("01_RawData/Breast", pattern = "cct$|iTRAQ", full.names = TRUE))
  f <- list.files("01_RawData/Breast", pattern = "iTRAQ.*cct$", full.names = TRUE)[1]
  if(!is.na(f)){ h <- strsplit(readLines(f, n = 1), "\t")[[1]]; cat("  samples in raw .cct:", length(h) - 1, "\n") }
})

# 4. WHICH DMP / DEP FILES EXIST, AND ARE THEY THE SAME? --------------------------
safe("4. DMP / DEP file locations", {
  fs <- c("03_Results/Epigenomics/brca_DMPs.csv","03_DifferentialAnalysis/Breast/brca_DMPs.csv",
          "results/step2_2_methylation_diff/breast_methylation_DMPs_filtered.csv",
          "03_Results/Proteomics/brca_DEPs.csv","03_DifferentialExpression/Breast/brca_DEPs.csv",
          "03_Results/Transcriptomics/brca_DEGs.csv","03_Results/Epigenomics/prad_DMPs.csv","03_Results/Transcriptomics/prad_DEGs.csv")
  for(f in fs){ info(f); if(file.exists(f)) cat("     rows:", length(readLines(f)) - 1, " md5:", tools::md5sum(f), "\n") }
})

# 5. RNA/meth saved objects & the brca_meth.rds save ------------------------------
safe("5. Key raw / processed files (timestamps)", {
  for(f in c("01_RawData/Breast/clinical_brca.csv","01_RawData/Breast/brca_rna.rds","01_RawData/Breast/brca_meth.rds",
             "01_RawData/Breast/brca_protein_raw.csv","01_RawData/Breast/brca_protein_early.csv",
             "02_Normalized/Breast/brca_protein_normalized.csv","02_ProcessedData/Breast/brca_protein_normalized.csv",
             "02_ProcessedData/Breast/brca_meth_normalized.csv","02_ProcessedData/Breast/dds_brca.rds")) info(f)
})
safe("5b. packages / versions", {
  for(p in c("TCGAbiolinks","sesameData","sesame","SummarizedExperiment","DESeq2","limma","minfi","MOFA2","mixOmics"))
    cat(sprintf("  %-22s %s\n", p, tryCatch(as.character(packageVersion(p)), error = function(e) "not installed")))
})

# 6. SEARCH YOUR SAVED R HISTORY FOR THE MISSING COMMANDS -------------------------
safe("6. .Rhistory search (RStudio: Tools > Global Options > General > may keep history)", {
  hf <- unique(c(".Rhistory", file.path(Sys.getenv("USERPROFILE"), "Documents", ".Rhistory"), file.path(Sys.getenv("USERPROFILE"), ".Rhistory"),
                 "C:/Users/USER/Documents/MyThesis_MultiOmics/.Rhistory"))
  hf <- hf[file.exists(hf)]; cat("  history files found:", hf, "\n")
  pat <- "sesameData|saveRDS\\(brca_meth|brca_protein_raw|value *=|sample_barcode|pivot_longer|melt\\(|PRIDE|brca_DMPs|brca_DEPs"
  for(h in hf){ l <- readLines(h, warn = FALSE); hit <- grep(pat, l); cat("  ", h, ":", length(hit), "matching lines\n"); print(l[hit]) }
})
safe("6b. Scripts saved anywhere in the project mentioning the missing writers", {
  fs <- list.files(".", pattern = "\\.(R|Rmd|txt)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  for(f in fs){ l <- tryCatch(readLines(f, warn = FALSE), error = function(e) character()); 
    hit <- grep("sample_barcode|brca_protein_raw.csv|saveRDS\\(brca_meth|sesameData::|BiocManager::install\\(.sesame", l)
    if(length(hit)) cat("  ", f, "-> lines", head(hit, 8), "\n") }
})
sink()
cat("\nDONE. Send 00_Scripts/confirm_outstanding_output.txt\n")
