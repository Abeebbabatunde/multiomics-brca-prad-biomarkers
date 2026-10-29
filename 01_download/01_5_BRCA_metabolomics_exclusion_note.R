# ================================================
# OPTION A: DOCUMENT BREAST METABOLOMICS
# EXCLUSION AND MOVE TO STEP 1.6
# ================================================

cat("========================================\n")
cat("   DOCUMENTING DECISION\n")
cat("========================================\n\n")

# -----------------------------------------------
# Create a note file documenting the decision
# -----------------------------------------------
exclusion_note <- data.frame(
  Decision   = "Breast cancer metabolomics excluded",
  Reason     = paste(
    "No publicly available breast cancer",
    "metabolomics dataset was found in a",
    "compatible format during data collection.",
    "MTBLS1 was identified as a diabetes study",
    "not breast cancer. GEO datasets GSE109856",
    "and GSE81827 contained empty data matrices.",
    "LinkedOmics TCGA-BRCA metabolomics file",
    "was not accessible at time of download."
  ),
  Thesis_Statement = paste(
    "Breast cancer metabolomics data was not",
    "available in a publicly accessible and",
    "computationally compatible format.",
    "Therefore metabolomics analysis was",
    "conducted for prostate cancer only.",
    "This represents a study limitation",
    "and future work should incorporate",
    "breast cancer metabolomics when data",
    "becomes available."
  ),
  Date = Sys.Date(),
  stringsAsFactors = FALSE
)

write.csv(
  exclusion_note,
  "03_Results/Metabolomics/breast_metabolomics_exclusion_note.csv",
  row.names = FALSE
)

cat("✅ Exclusion note saved\n\n")

# -----------------------------------------------
# Create empty placeholder so startup script
# does not flag it as missing
# -----------------------------------------------
brca_met_excluded <- TRUE
saveRDS(
  brca_met_excluded,
  "01_RawData/Breast/brca_metabolomics_excluded.rds"
)

cat("✅ Placeholder saved\n\n")

# -----------------------------------------------
# Update your thesis methodology text
# -----------------------------------------------
cat("========================================\n")
cat("   YOUR THESIS METHODOLOGY TEXT\n")
cat("========================================\n\n")

cat("Copy this into your thesis:\n\n")
cat("'Metabolomics data for breast cancer\n")
cat(" could not be obtained in a publicly\n")
cat(" accessible and computationally\n")
cat(" compatible format despite systematic\n")
cat(" searches of MetaboLights (MTBLS1,\n")
cat(" MTBLS522), NCBI GEO (GSE109856,\n")
cat(" GSE81827), and LinkedOmics databases.\n")
cat(" Consequently, metabolomics analysis\n")
cat(" was performed exclusively for prostate\n")
cat(" cancer using data from MetaboLights\n")
cat(" MTBLS103. Breast cancer multi-omics\n")
cat(" integration was performed using three\n")
cat(" omics layers: epigenomics, trans-\n")
cat(" criptomics, and proteomics. This\n")
cat(" represents a study limitation and\n")
cat(" future studies should incorporate\n")
cat(" breast cancer metabolomics data.'\n\n")

# -----------------------------------------------
# Show updated progress
# -----------------------------------------------
cat("========================================\n")
cat("   UPDATED PROGRESS\n")
cat("========================================\n\n")

steps_status <- list(
  list("1.1","Breast clinical",
       "01_RawData/Breast/clinical_brca.csv"),
  list("1.2","Breast RNA-seq",
       "01_RawData/Breast/brca_rna.rds"),
  list("1.3","Breast methylation",
       "01_RawData/Breast/brca_meth.rds"),
  list("1.4","Breast proteomics",
       "01_RawData/Breast/brca_protein_early.csv"),
  list("1.5","Breast metabolomics",
       "01_RawData/Breast/brca_metabolomics_excluded.rds"),
  list("1.6","Prostate clinical",
       "01_RawData/Prostate/clinical_prad.csv"),
  list("1.7","Prostate RNA-seq",
       "01_RawData/Prostate/prad_rna.rds"),
  list("1.8","Prostate methylation",
       "01_RawData/Prostate/prad_meth.rds"),
  list("1.9","Prostate proteomics",
       "01_RawData/Prostate/prad_protein_early.csv"),
  list("1.10","Prostate metabolomics",
       "01_RawData/Prostate/prad_metabolites_raw.csv")
)

done_count <- 0
next_step  <- NULL

cat(sprintf("  %-6s %-30s %s\n",
            "STEP","NAME","STATUS"))
cat(paste(rep("-", 55), collapse=""), "\n")

for(s in steps_status){
  done   <- file.exists(s[[3]])
  status <- ifelse(done, "✅ DONE", "⬜ PENDING")
  cat(sprintf("  %-6s %-30s %s\n",
              s[[1]], s[[2]], status))
  if(done) done_count <- done_count + 1
  if(!done && is.null(next_step)) next_step <- s
}

cat(paste(rep("-", 55), collapse=""), "\n")
cat(sprintf("  Progress: %d / %d\n",
            done_count, length(steps_status)))
cat("========================================\n")

if(!is.null(next_step)){
  cat(sprintf("\n  ▶️  NEXT STEP: %s — %s\n",
              next_step[[1]],
              next_step[[2]]))
  cat("========================================\n")
}