# STEP 1.6 — Prostate clinical data (TCGA-PRAD). Early = pathologic T1/T2.
# Console: 188 patients; T2a 13, T2b 10, T2c 165 (unique T stages: T2a, T2b, T2c). Needs library(TCGAbiolinks). A "Timeout was reached [api.gdc.cancer.gov]" error occurred once; re-run.
cat("Downloading prostate clinical data...\n")
clinical_prad_full <- GDCquery_clinic("TCGA-PRAD", type = "clinical")
clinical_prad <- clinical_prad_full[
  clinical_prad_full$ajcc_pathologic_t %in% c("T1","T1a","T1b","T1c","T2","T2a","T2b","T2c"), ]
early_barcodes_prad <- clinical_prad$submitter_id
cat("Early prostate patients:", length(early_barcodes_prad), "\n")
print(table(clinical_prad$ajcc_pathologic_t))
write.csv(clinical_prad, "01_RawData/Prostate/clinical_prad.csv", row.names = FALSE)
