# Step 0.5 add-on: are BRCA patients with proteomics representative of the BRCA manifest?
# Reconstructed from the code block that was run (output: chi-square / Fisher / tableone, shown in the console).
# NOTE: the header comment says "and the 31 with all 3 omics", but the code compares has_proteomics (76 vs 727) only.
library(dplyr); library(readr); library(tableone)
manifest <- read_csv("00.5_clinical_metadata_harmonization/sample_manifest_harmonized.csv", show_col_types = FALSE) %>%
  filter(project == "BRCA")
table(manifest$has_proteomics); table(manifest$eligible_for_integration)
print(table(manifest$ajcc_pathologic_stage, manifest$has_proteomics))
chisq.test(table(manifest$ajcc_pathologic_stage, manifest$has_proteomics))
fisher.test(table(manifest$ajcc_pathologic_stage, manifest$has_proteomics), simulate.p.value = TRUE, B = 10000)  # unseeded in the run
print(table(manifest$tumor_grade, manifest$has_proteomics))
print(table(manifest$vital_status, manifest$has_proteomics))
chisq.test(table(manifest$vital_status, manifest$has_proteomics))
tab1 <- CreateTableOne(vars = c("ajcc_pathologic_stage","tumor_grade","vital_status"), strata = "has_proteomics", data = manifest)
print(tab1, showAllLevels = TRUE)
