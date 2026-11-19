# 00b_audit_pairing.R
# Step 12.0b (Step A): audit of patient pairing in the cohort manifests (read-only; changes no data).
# Checks the rebuild manifests (12_rebuild_inputs/cohort_*_rebuild.csv) for BRCA and PRAD:
#  - sample-type codes in the barcodes (tumour = 01, normal = 11) and patient ID = first 12 characters
#  - how many tumour and normal vials each patient has in the processed matrices (the cohort script picks the FIRST sorted barcode)
#  - whether the same patients and the same barcodes are used in every layer
# Edit project_dir only.

project_dir <- "C:/Users/USER/Documents/MyThesis_MultiOmics"
setwd(project_dir)
fix  <- function(x) gsub(".", "-", x, fixed = TRUE)
code <- function(b) substr(fix(b), 14, 15); pat <- function(b) substr(fix(b), 1, 12)

check_manifest <- function(file, layers, matrices) {
  man <- read.csv(file, stringsAsFactors = FALSE)
  cat("\n====", file, "====\n")
  cat("Rows:", nrow(man), " patients:", length(unique(man$patient)), "\n")
  print(table(man$group))
  stopifnot(all(table(man$patient) == 2))                                   # one Normal and one Early row per patient
  for (l in layers) {
    b <- fix(man[[l]])
    cat(l, ": patient ID matches barcode:", all(pat(b) == man$patient),
        "| codes (Normal):", paste(unique(code(b[man$group == "Normal"])), collapse = ","),
        "| codes (Early):", paste(unique(code(b[man$group == "Early"])), collapse = ","), "\n")
  }
  for (l in names(matrices)) {                                              # vials per patient and class in the full matrix
    h <- fix(matrices[[l]]); h <- h[grepl("^TCGA", h)]
    d <- data.frame(patient = pat(h), cls = ifelse(code(h) == "11", "Normal", ifelse(code(h) == "01", "Tumour", "Other")))
    d <- d[d$patient %in% man$patient & d$cls != "Other", ]
    tab <- table(d$patient, d$cls)
    cat(l, ": patients with >1 tumour vial:", sum(tab[, "Tumour"] > 1), "; with >1 normal vial:", sum(tab[, "Normal"] > 1), "\n")
  }
}
library(data.table)
hdr <- function(f) names(fread(f, nrows = 0))
check_manifest("12_rebuild_inputs/cohort_BRCA_rebuild.csv", c("meth", "rna"),
               list(methylation = hdr("02_ProcessedData/Breast/brca_meth_normalized.csv"),
                    rna = as.data.frame(readRDS("02_ProcessedData/Breast/brca_rna_meta.rds"))$barcode))
check_manifest("12_rebuild_inputs/cohort_PRAD_rebuild.csv", c("meth", "rna"),
               list(methylation = hdr("02_ProcessedData/Prostate/prad_meth_normalized.csv"),
                    rna = as.data.frame(readRDS("02_ProcessedData/Prostate/prad_rna_meta.rds"))$barcode))
cat("\nIf any patient has more than one tumour or normal vial, record that the first sorted barcode was used.\n")
