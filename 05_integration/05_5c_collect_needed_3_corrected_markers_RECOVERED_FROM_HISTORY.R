## 05_5c  RECOVERED FROM RStudio HISTORY -- collect_needed_3.R, run 2026-10-19 08:08:00-08:08:30
## Wrote PRAD_genes_in_all_layers_corrected.csv and classic_marker_validation_exact_corrected.csv (corrected DMP table,
## whole-symbol matching) plus needed_outputs_3.txt. Section 1 reads the stage of the five unmodelled BRCA patients.
## The 00_STARTUP.R creation code that precedes it in the same history window (08:07:40) is unrelated and left out.
## Lines are exactly the history lines; RStudio history does not store indentation or blank lines.
# ============================================================================
# collect_needed_3.R
# Recomputes the PRAD overlap, DMP share and classical-marker results on the
# CORRECTED DMP table (79,679 DMPs), and checks the stage of the five BRCA
# patients that were paired in RNA and methylation but not modelled.
#
# It does NOT rerun DIABLO or MOFA and does not overwrite any existing file.
# New files written to the project folder:
#   PRAD_genes_in_all_layers_corrected.csv
#   classic_marker_validation_exact_corrected.csv
# Run:  source("collect_needed_3.R")   then upload needed_outputs_3.txt
# ============================================================================
setwd("C:/Users/USER/Documents/MyThesis_MultiOmics")
options(width = 170, max.print = 300, warn = 1)
suppressPackageStartupMessages({
library(data.table); library(org.Hs.eg.db)
library(minfi); library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
})
b <- "05_FinalResults/Step5.0_Integration"
while (sink.number() > 0) sink()
sink("needed_outputs_3.txt", split = TRUE)
sec <- function(title, expr) {
cat("\n\n############################################################\n")
cat("##", title, "\n")
cat("############################################################\n")
tryCatch(expr, error = function(e) cat("SECTION FAILED:", conditionMessage(e), "\n"))
invisible(NULL)
}
sec("1. STAGE OF THE FIVE BRCA PATIENTS PAIRED IN RNA AND METHYLATION BUT NOT MODELLED", {
m <- fread("00.5_clinical_metadata_harmonization/sample_manifest_harmonized.csv")
five <- c("TCGA-A7-A13E", "TCGA-E9-A1N6", "TCGA-E9-A1NF", "TCGA-BH-A0C3", "TCGA-A7-A13G")
cat("present in the harmonised manifest:", sum(five %in% m$patient_id), "of 5\n")
print(m[patient_id %in% five, .(patient_id, ajcc_pathologic_stage, vital_status, has_epigenomics, has_transcriptomics, has_proteomics)])
cat("\nBRCA stages in the harmonised manifest (803 patients):\n")
print(table(m[project == "BRCA", ajcc_pathologic_stage], useNA = "ifany"))
})
sec("2. PRAD: OVERLAP, DMP SHARE AND MARKERS ON THE CORRECTED DMP TABLE", {
d  <- fread("03_Results/Epigenomics/prad_DMPs_corrected_results.csv")
dm <- d[adj.P.Val < 0.05 & !is.na(delta_beta_corrected) & abs(delta_beta_corrected) >= 0.1]
cat("corrected DMPs:", nrow(dm), "(1.0 reports 79,679) | higher in cancer:", sum(dm$delta_beta_corrected > 0),
"| lower:", sum(dm$delta_beta_corrected < 0), "\n")
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
g   <- ann[match(dm$probe_id, rownames(ann)), "UCSC_RefGene_Name"]
cat("DMP probes found in the annotation:", sum(dm$probe_id %in% rownames(ann)), "of", nrow(dm), "\n")
meth_set <- unique(na.omit(trimws(unlist(strsplit(as.character(g), ";")))))
deg <- fread("03_Results/Transcriptomics/prad_DEGs.csv")
rna_set <- unique(na.omit(mapIds(org.Hs.eg.db, sub("\\..*", "", deg$ensembl_gene_id), "SYMBOL", "ENSEMBL", multiVals = "first")))
common <- intersect(meth_set, rna_set)
cat("\nmethylation genes:", length(meth_set), "| transcriptome genes:", length(rna_set), "| overlap:", length(common), "\n")
cat("(original table: 11,490 methylation genes, 3,776 transcriptome genes, overlap 1,798)\n")
fwrite(data.frame(gene = common), "PRAD_genes_in_all_layers_corrected.csv")
sel <- fread(file.path(b, "PRAD/meth_rna_17pairs/DIABLO_selected_features.csv"))
mp  <- unique(sel$feature[sel$block == "meth"])
cat("\nselected PRAD methylation probes also in the corrected DMP set:", sum(mp %in% dm$probe_id), "of", length(mp),
"(original table: 10 of 20)\n")
cat("\nclassical markers, whole-symbol matching, corrected DMP set:\n")
sym <- strsplit(as.character(g), ";")
out <- rbindlist(lapply(c("GSTP1", "RASSF1", "APC", "PTGS2"), function(mk) {
h <- vapply(sym, function(x) mk %in% x, logical(1))
x <- dm[h, .(probe_id, delta_beta = delta_beta_corrected, logFC_Mvalue = logFC, adj.P.Val)]
x[, `:=`(gene = g[h], marker = mk)]
x[order(adj.P.Val)]
}))
fwrite(out, "classic_marker_validation_exact_corrected.csv")
print(out[, .(CpGs = .N, smallest_adjP = min(adj.P.Val), higher_in_cancer = sum(delta_beta > 0),
lower_in_cancer = sum(delta_beta < 0), min_delta = min(delta_beta), max_delta = max(delta_beta)), by = marker])
cat("(original table: GSTP1 7, RASSF1 17, APC 17, PTGS2 7 CpGs)\n")
})
cat("\n\nDONE. Upload needed_outputs_3.txt.\n")
while (sink.number() > 0) sink()
cat("DONE. Output written to", file.path(getwd(), "needed_outputs_3.txt"), "\n")
