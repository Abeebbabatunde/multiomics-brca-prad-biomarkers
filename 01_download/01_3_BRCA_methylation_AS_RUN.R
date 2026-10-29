# STEP 1.3 — Breast methylation (450K beta) for the early-stage patient barcodes (patient barcodes -> all sample types).
# Console: "Merging 649 files" (matches the 649-sample metadata seen in Step 3.2 diagnostics).
# Failure seen once at GDCprepare: "sesameData package is needed for this function to work. Please install it."
#   -> install sesameData (the install command was not pasted) and re-run GDCprepare.
# Another time GDCquery failed ("GDC server down ... SSL connect error"): re-run later.
# Integrity checks pasted: GDCcheck(query_brca_meth, directory="01_RawData/Breast/").
options(timeout = 7200)
if(!exists("early_barcodes_brca")){
  clinical_brca       <- read.csv("01_RawData/Breast/clinical_brca.csv")
  early_barcodes_brca <- clinical_brca$submitter_id
}
query_brca_meth <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "DNA Methylation",
  data.type     = "Methylation Beta Value",
  platform      = "Illumina Human Methylation 450",
  barcode       = early_barcodes_brca
)
GDCdownload(query_brca_meth, method = "api", files.per.chunk = 4, directory = "01_RawData/Breast/")
brca_meth        <- GDCprepare(query_brca_meth, directory = "01_RawData/Breast/")
meth_matrix_brca <- assay(brca_meth)
meth_meta_brca   <- colData(brca_meth)
cat("✅ Probes :", nrow(meth_matrix_brca), "\n"); cat("✅ Samples:", ncol(meth_matrix_brca), "\n")
# NOT pasted: the saveRDS(brca_meth, "01_RawData/Breast/brca_meth.rds") line (the file is read by later steps; the PRAD analogue is saveRDS(prad_meth, ...)).
