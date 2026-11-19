# Sex and gender check for the two rebuild cohorts (Step 12.0)
# Run from the project root (the folder that holds 01_RawData and 02_ProcessedData),
# after 00_define_cohorts_and_inputs.R has written the two cohort files.
# Set out_dir to the folder where cohort_BRCA_rebuild.csv and cohort_PRAD_rebuild.csv were saved.
out_dir <- "12_rebuild_inputs"

man_b <- read.csv(file.path(out_dir, "cohort_BRCA_rebuild.csv"), stringsAsFactors = FALSE)
man_p <- read.csv(file.path(out_dir, "cohort_PRAD_rebuild.csv"), stringsAsFactors = FALSE)
rm_b  <- as.data.frame(readRDS("02_ProcessedData/Breast/brca_rna_meta.rds"))
rm_p  <- as.data.frame(readRDS("02_ProcessedData/Prostate/prad_rna_meta.rds"))
clin_b <- read.csv("01_RawData/Breast/clinical_brca.csv", stringsAsFactors = FALSE)

report <- function(label, patients, meta, cols) {
  cat("\n==", label, "==\n")
  cat("Patients in the cohort:", length(unique(patients)), "\n")
  for (cl in intersect(cols, names(meta))) {
    d <- meta[!duplicated(meta$patient) & meta$patient %in% patients, , drop = FALSE]
    cat("Column '", cl, "' (", nrow(d), " cohort patients found in this table):\n", sep = "")
    print(table(d[[cl]], useNA = "ifany"))
  }
}

report("Breast, RNA metadata",    man_b$patient, rm_b, c("gender", "sex_at_birth"))
report("Prostate, RNA metadata",  man_p$patient, rm_p, c("gender", "sex_at_birth"))

# The clinical file may name its patient column differently
pc <- intersect(c("patient", "bcr_patient_barcode", "submitter_id"), names(clin_b))[1]
clin_b$patient <- clin_b[[pc]]
report("Breast, clinical file",   man_b$patient, clin_b, c("gender", "sex_at_birth"))
