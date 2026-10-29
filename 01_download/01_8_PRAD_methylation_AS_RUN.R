# STEP 1.8 — Prostate methylation (450K beta values) for the 188 early patients.
# Check afterwards: rows > 400000 and columns > 50.
options(timeout = 7200)
if (!exists("early_barcodes_prad")) {
  clinical_prad <- read.csv("01_RawData/Prostate/clinical_prad.csv")
  early_barcodes_prad <- clinical_prad$submitter_id
}
query_prad_meth <- GDCquery(project = "TCGA-PRAD", data.category = "DNA Methylation",
  data.type = "Methylation Beta Value", platform = "Illumina Human Methylation 450", barcode = early_barcodes_prad)
GDCdownload(query_prad_meth, method = "api", files.per.chunk = 4, directory = "01_RawData/Prostate/")
prad_meth        <- GDCprepare(query_prad_meth, directory = "01_RawData/Prostate/")
meth_matrix_prad <- assay(prad_meth)
meth_meta_prad   <- colData(prad_meth)
saveRDS(prad_meth, "01_RawData/Prostate/prad_meth.rds")
write.csv(meth_matrix_prad, "01_RawData/Prostate/prad_meth_matrix.csv")
saveRDS(meth_meta_prad, "01_RawData/Prostate/prad_meth_meta.rds")
