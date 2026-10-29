## 00_5  Paired Normal-vs-Early BRCA manifest + hand-entered attrition table
## Recovered verbatim from the RStudio history, run 2026-09-20 16:32:48. It writes the two files that Step 7 reads
## (sample_manifest_paired_normal_early.csv) and the attrition_report_paired_design.csv.
## Needs 05_0a to have written cohort_3omics.csv (16 pairs) and cohort_meth_rna.csv (51 pairs) first.
## NOTE: the numbers in `n` are typed by hand in the original (80, 64, 28, 16, 51). 16 and 51 match the cohorts;
## 80/64/28 should match the sapply(P, length) printout in cohort_check.txt (not yet compared).
setwd("C:/Users/USER/Documents/MyThesis_MultiOmics"); library(data.table)
b  <- "05_FinalResults/Step5.0_Integration/BRCA"; d <- "00.5_clinical_metadata_harmonization"
mf <- fread(file.path(d, "sample_manifest_harmonized.csv"))
c3 <- fread(file.path(b, "cohort_3omics.csv"));   c3[, cohort := "3omics_16pairs"]
c2 <- fread(file.path(b, "cohort_meth_rna.csv")); c2[, cohort := "meth_rna_51pairs"]
pm <- merge(rbind(c3, c2, fill = TRUE),
            mf[, .(patient_id, ajcc_pathologic_stage, vital_status, days_to_death, days_to_last_follow_up)],
            by.x = "patient", by.y = "patient_id", all.x = TRUE)
fwrite(pm, file.path(d, "sample_manifest_paired_normal_early.csv"))
att <- data.frame(
  step = c("BRCA patients with Normal + Stage I-II tumor: RNA", "...: methylation", "...: proteomics (30 pairs file)",
           "All 3 omics (final exploratory cohort)", "Methylation + RNA (final main cohort)"),
  n = c(80, 64, 28, 16, 51))
fwrite(att, file.path(d, "attrition_report_paired_design.csv")); print(att)
