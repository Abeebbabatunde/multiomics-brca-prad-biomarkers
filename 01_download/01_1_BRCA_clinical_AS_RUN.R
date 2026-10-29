# STEP 1.1 — Breast clinical data; Early = AJCC stage I/IA/IB/II/IIA/IIB (804 patients per later Step 5 tables).
# Assembled from pasted fragments (original run order not shown). Needs library(TCGAbiolinks).
library(TCGAbiolinks)
cat("Downloading breast cancer clinical data...\n")
clinical_brca_full <- GDCquery_clinic("TCGA-BRCA", type = "clinical")

clinical_brca <- clinical_brca_full[
  clinical_brca_full$ajcc_pathologic_stage %in% c(
    "Stage I","Stage IA","Stage IB",
    "Stage II","Stage IIA","Stage IIB"
  ), ]
early_barcodes_brca <- clinical_brca$submitter_id
cat("✅ Early-stage breast cancer patients:", length(early_barcodes_brca), "\n")
cat("Stage distribution:\n"); print(table(clinical_brca$ajcc_pathologic_stage))

# write.csv() first failed: "unimplemented type 'list'" — column 31 `sites_of_involvement` is a list column.
# Fix pasted (flatten list column), then write:
clinical_brca$sites_of_involvement <- sapply(
  clinical_brca$sites_of_involvement,
  function(x) paste(unlist(x), collapse = "; ")
)
write.csv(clinical_brca, "01_RawData/Breast/clinical_brca.csv", row.names = FALSE)
cat("✅ SAVED: 01_RawData/Breast/clinical_brca.csv\n")
# (An alternative dplyr version, mutate(across(where(is.list), ~ sapply(., paste, collapse="; "))), was also pasted.)
# An earlier integrity-check script refers to "01_RawData/Breast/clinical_early_stage_brca.csv" (earlier file name); see archive/.
