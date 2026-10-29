# STEP 1.7 — Prostate RNA-seq (STAR counts) for the 188 early patients. Needs TCGAbiolinks, SummarizedExperiment.
# (assay() failed once with "could not find function" until SummarizedExperiment was loaded.)
options(timeout = 7200)
if (!exists("early_barcodes_prad")) {
  clinical_prad <- read.csv("01_RawData/Prostate/clinical_prad.csv")
  early_barcodes_prad <- clinical_prad$submitter_id
}
query_prad_rna <- GDCquery(project = "TCGA-PRAD", data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts", barcode = early_barcodes_prad)
GDCdownload(query_prad_rna, method = "api", files.per.chunk = 6, directory = "01_RawData/Prostate/")
prad_rna        <- GDCprepare(query_prad_rna, directory = "01_RawData/Prostate/")
rna_counts_prad <- assay(prad_rna, "unstranded")
rna_meta_prad   <- colData(prad_rna)
saveRDS(prad_rna, "01_RawData/Prostate/prad_rna.rds")
write.csv(rna_counts_prad, "01_RawData/Prostate/prad_rna_counts.csv")
# (saveRDS(rna_meta_prad) failed once with "object not found" because the object was missing in that session; re-run the colData() line first.)
