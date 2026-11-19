# 03_compare_old_and_rebuilt_panels.R
# Step 12.3 (Objective 2): compares the locked thesis panels with the rebuilt panels and the selection frequencies.
# Reads the Step 10 panel files (same paths as 11_pre_submission_reanalyses/01_nested_patient_cv.R) and the files from 12.1.
# Edit project_dir only.

project_dir <- "C:/Users/USER/Documents/MyThesis_MultiOmics"
setwd(project_dir)
cfg <- list(
  PRAD_methylation = list(panel = "10_FinalResults/Step10_FinalPanel/Step10_2_PRAD_final_panel.csv", omics = "methylation", ensembl = FALSE, x = "12_rebuild_inputs/X_methylation_PRAD_rebuild.rds"),
  BRCA_methylation = list(panel = "10_FinalResults/Step10_FinalPanel/Step10_1_BRCA_final_panel.csv", omics = "methylation", ensembl = FALSE, x = "12_rebuild_inputs/X_methylation_BRCA_rebuild.rds"),
  BRCA_RNA         = list(panel = "10_FinalResults/Step10_FinalPanel/Step10_1_BRCA_final_panel.csv", omics = "RNA", ensembl = TRUE, x = "12_rebuild_inputs/X_rna_BRCA_rebuild.rds"),
  PRAD_RNA         = list(panel = "10_FinalResults/Step10_FinalPanel/Step10_2_PRAD_final_panel.csv", omics = "RNA", ensembl = TRUE, x = "12_rebuild_inputs/X_rna_PRAD_rebuild.rds"))

out <- do.call(rbind, lapply(names(cfg), function(nm) {
  c1 <- cfg[[nm]]
  pn  <- read.csv(c1$panel, stringsAsFactors = FALSE)
  old <- unique(as.character(pn$feature[pn$omics == c1$omics]))
  fd  <- read.csv(paste0(nm, "_selection_frequency_rebuild.csv"), stringsAsFactors = FALSE)
  new <- read.csv(paste0(nm, "_rebuilt_panel.csv"), stringsAsFactors = FALSE)
  id  <- if (c1$ensembl) fd$symbol else fd$feature
  f_old <- fd$selection_frequency[match(old, id)]; f_old[is.na(f_old)] <- 0       # 0 = never selected
  cn <- colnames(readRDS(c1$x))                                   # features that were actually available to the rebuild
  tested_id <- if (c1$ensembl) as.character(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = sub("\\.[0-9]+$", "", cn),
                 keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")) else cn
  n_tested <- sum(old %in% tested_id)
  newid <- if (c1$ensembl) new$symbol else new$feature
  data.frame(layer = nm, old_panel_size = length(old), rebuilt_panel_size = nrow(new),
             in_both = length(intersect(old, newid)),
             old_features_median_frequency = median(f_old),
             old_features_tested = n_tested, old_features_not_tested = length(old) - n_tested,
              old_features_tested_never_selected = sum(f_old == 0 & old %in% tested_id))
}))
print(out, row.names = FALSE)
write.csv(out, "rebuild_03_old_vs_rebuilt.csv", row.names = FALSE)
cat("\nA rebuilt panel that overlaps little with the locked panel supports the thesis finding that the locked panel was one of several possible solutions.\n")
