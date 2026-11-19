# 00_check_patient_groups.R
# Step 12.0: checks that the rebuild inputs are ready and that every patient has BOTH a tumour and a normal sample,
# so that patient-grouped folds are possible. It changes no data. Run after 00_define_cohorts_and_inputs.R.
# Edit project_dir only.

project_dir <- "C:/Users/USER/Documents/MyThesis_MultiOmics"
setwd(project_dir)

layers <- list(
  PRAD_methylation = list(X = "12_rebuild_inputs/X_methylation_PRAD_rebuild.rds", meta = "12_rebuild_inputs/meta_PRAD_methylation_rebuild.csv", transpose = FALSE),
  BRCA_methylation = list(X = "12_rebuild_inputs/X_methylation_BRCA_rebuild.rds", meta = "12_rebuild_inputs/meta_BRCA_methylation_rebuild.csv", transpose = FALSE),
  BRCA_RNA         = list(X = "12_rebuild_inputs/X_rna_BRCA_rebuild.rds",         meta = "12_rebuild_inputs/meta_BRCA_rna_rebuild.csv", transpose = FALSE),
  PRAD_RNA         = list(X = "12_rebuild_inputs/X_rna_PRAD_rebuild.rds",         meta = "12_rebuild_inputs/meta_PRAD_rna_rebuild.csv", transpose = FALSE))

report <- do.call(rbind, lapply(names(layers), function(nm) {
  cfg <- layers[[nm]]
  stopifnot(file.exists(cfg$X), file.exists(cfg$meta))
  X <- readRDS(cfg$X); if (cfg$transpose) X <- t(X)
  meta <- read.csv(cfg$meta, stringsAsFactors = FALSE)
  meta <- meta[match(rownames(X), meta$sample), ]
  tab <- table(meta$patient, meta$class)
  data.frame(layer = nm, samples = nrow(X), features = ncol(X), patients = nrow(tab),
             tumour = sum(meta$class == "Tumour"), normal = sum(meta$class == "Normal"),
             patients_with_both = sum(tab[, "Normal"] > 0 & tab[, "Tumour"] > 0),
             patients_with_more_than_one_sample_per_class = sum(tab > 1),
             missing_values = sum(is.na(X)))
}))
print(report, row.names = FALSE)
if (any(report$patients != report$patients_with_both)) warning("Some patients lack a tumour or a normal sample.")
cat("\nCheck the numbers against the thesis: the pair counts printed by 00_define_cohorts_and_inputs.R (thesis: breast 51, prostate 17; the matrices hold two samples per patient).\n")
write.csv(report, "12_rebuild_inputs/rebuild_00_input_check.csv", row.names = FALSE)
