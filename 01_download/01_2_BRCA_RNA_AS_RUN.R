# STEP 1.2 — Breast RNA-seq (STAR - Counts) for the early-stage patient barcodes.
# NB: `barcode = early_barcodes_brca` are PATIENT barcodes, so ALL sample types of those patients are retrieved
# (tumour, solid tissue normal, and any others). This explains the sample-type composition question in Step 2.1.
# Assembled from pasted fragments. Needs TCGAbiolinks, SummarizedExperiment.
options(timeout = 7200)
if(!exists("early_barcodes_brca")){
  clinical_brca       <- read.csv("01_RawData/Breast/clinical_brca.csv")
  early_barcodes_brca <- clinical_brca$submitter_id
  cat("✅ Barcodes reloaded:", length(early_barcodes_brca), "\n")
}
cat("Building download query...\n")
query_brca_rna <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode       = early_barcodes_brca
)
total <- nrow(getResults(query_brca_rna)); cat("Files to download:", total, "\n")
GDCdownload(query = query_brca_rna, method = "api", files.per.chunk = 6, directory = "01_RawData/Breast/")
brca_rna        <- GDCprepare(query_brca_rna, directory = "01_RawData/Breast/")
rna_counts_brca <- assay(brca_rna, "unstranded")
rna_meta_brca   <- colData(brca_rna)
saveRDS(brca_rna, "01_RawData/Breast/brca_rna.rds")
write.csv(rna_counts_brca, "01_RawData/Breast/brca_rna_counts.csv")
saveRDS(rna_meta_brca, "01_RawData/Breast/brca_rna_metadata.rds")
