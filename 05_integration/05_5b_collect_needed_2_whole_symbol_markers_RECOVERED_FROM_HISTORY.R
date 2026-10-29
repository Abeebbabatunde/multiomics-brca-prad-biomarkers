## 05_5b  RECOVERED FROM RStudio HISTORY -- collect_needed_2.R, run 2026-10-19 02:27:06-02:28:00
## Section 6 of this script is the whole-symbol classical-marker script: it wrote classic_marker_validation_exact.csv
## (original prad_DMPs.csv). Sections 1-5 are read-only checks (attrition files, BRCA methylation pathway CSVs,
## PRAD DMP counts, BRCA 56 vs 51 pairs, PRAD RNA scale) and wrote only needed_outputs_2.txt.
## Lines are exactly the history lines; RStudio history does not store indentation or blank lines.
# ============================================================================
# collect_needed_2.R
# Second, smaller collection for the items that are still open.
# Run:  source("collect_needed_2.R")   then upload needed_outputs_2.txt.
# Nothing is modified except one NEW file written by section 6
# (classic_marker_validation_exact.csv, in the project folder).
# ============================================================================
setwd("C:/Users/USER/Documents/MyThesis_MultiOmics")
options(width = 170, max.print = 300, warn = 1)
suppressPackageStartupMessages(library(data.table))
b <- "05_FinalResults/Step5.0_Integration"
while (sink.number() > 0) sink()
sink("needed_outputs_2.txt", split = TRUE)
sec <- function(title, expr) {
cat("\n\n############################################################\n")
cat("##", title, "\n")
cat("############################################################\n")
tryCatch(expr, error = function(e) cat("SECTION FAILED:", conditionMessage(e), "\n"))
invisible(NULL)
}
ex   <- function(f) { ok <- file.exists(f); if (!ok) cat("MISSING:", f, "\n"); ok }
fix  <- function(x) gsub(".", "-", x, fixed = TRUE)
pat  <- function(x) substr(fix(x), 1, 12)
code <- function(x) substr(fix(x), 14, 15)
sec("1. THE PAIRED-DESIGN ATTRITION FILES (00.5 folder)", {
d <- "00.5_clinical_metadata_harmonization"
for (f in c("attrition_report.csv", "attrition_report_paired_design.csv")) {
cat("\n==", f, "==\n"); if (ex(file.path(d, f))) print(fread(file.path(d, f)))
}
for (f in c("sample_manifest_paired_normal_early.csv", "sample_manifest_harmonized.csv")) {
cat("\n==", f, "==\n")
if (!ex(file.path(d, f))) next
m <- fread(file.path(d, f)); cat(nrow(m), "rows;", ncol(m), "columns\n"); print(names(m)); print(head(m, 6))
for (cc in names(m)) if (is.character(m[[cc]]) && length(unique(m[[cc]])) <= 12) { cat("\n", cc, ":\n"); print(table(m[[cc]], useNA = "ifany")) }
}
})
sec("2. BRCA METHYLATION PATHWAY CSVs: how many significant terms", {
for (f in list.files("04_Pathway_Analysis/Epigenomics", pattern = "\\.csv$", full.names = TRUE)) {
d <- fread(f); cat("\n==", f, "==\n"); print(names(d))
pc <- intersect(c("padj", "p.adjust", "FDR"), names(d))[1]
cat("terms:", nrow(d), "| P exactly 0:", sum(d$pvalue == 0, na.rm = TRUE),
"| adjusted P < 0.05:", sum(d[[pc]] < 0.05, na.rm = TRUE), "| < 0.1:", sum(d[[pc]] < 0.1, na.rm = TRUE), "\n")
print(head(d[order(d[[pc]], d$pvalue), intersect(c("ID", "Description", "Count", "Size", "pvalue", pc), names(d)), with = FALSE], 12))
}
})
sec("3. PRAD DMP COUNTS WITH THE CORRECTED DELTA-BETA COLUMNS", {
d <- fread("03_Results/Epigenomics/prad_DMPs_corrected_results.csv")
cat("thesis text (1.0): 79,679 DMPs; 53,287 hypermethylated; 26,392 hypomethylated\n")
for (e in intersect(c("delta_beta", "recalc_delta_beta", "delta_beta_corrected"), names(d))) {
s <- d$adj.P.Val < 0.05 & abs(d[[e]]) >= 0.1
cat(sprintf("%-22s DMPs %d | higher in cancer %d | lower in cancer %d | missing values %d\n",
e, sum(s, na.rm = TRUE), sum(s & d[[e]] > 0, na.rm = TRUE), sum(s & d[[e]] < 0, na.rm = TRUE), sum(is.na(d[[e]]))))
}
cat("\nhow the three delta-beta columns differ:\n")
print(summary(d[, c("delta_beta", "recalc_delta_beta", "delta_beta_corrected"), with = FALSE]))
cat("probes where delta_beta and delta_beta_corrected differ by more than 0.001:",
sum(abs(d$delta_beta - d$delta_beta_corrected) > 0.001, na.rm = TRUE), "\n")
})
sec("4. BRCA: THE 56 PATIENTS PAIRED IN RNA AND METHYLATION VERSUS THE 51 MODELLED", {
paired <- function(ids) { p <- pat(ids); cd <- code(ids); intersect(unique(p[cd == "01"]), unique(p[cd == "11"])) }
rna  <- colnames(readRDS("01_RawData/Breast/brca_rna.rds"))
meth <- colnames(readRDS("02_ProcessedData/Breast/brca_meth_filtered_matrix.rds"))
both2 <- intersect(paired(rna), paired(meth))
e <- new.env(); load(file.path(b, "BRCA/prep_objects.RData"), envir = e)
p51 <- unique(e$A2$patient); p16 <- unique(e$A3$patient)
extra <- setdiff(both2, p51)
cat("paired in both layers:", length(both2), "| modelled:", length(p51), "| not modelled:", length(extra), "\n")
for (p in extra) {
cat("\n", p, "\n  RNA samples :", paste(rna[pat(rna) == p], collapse = ", "),
"\n  meth samples:", paste(meth[pat(meth) == p], collapse = ", "), "\n")
}
cat("\nthe 5 not modelled are in the eligible 31:",
sum(extra %in% readRDS("05_FinalResults/Step5_BRCA_3omics_patient_ids.rds")), "\n")
pf <- "01_RawData/Breast/brca_protein_sample_pairing.csv"
if (ex(pf)) {
pp <- fread(pf); cat("\nRPPA pairing file:", nrow(pp), "rows\n"); print(names(pp)); print(head(pp, 4))
pcol <- intersect(c("PatientID", "patient", "patient_id", "Patient"), names(pp))[1]
if (!is.na(pcol)) {
rp <- unique(pp[[pcol]])
cat("RPPA patients:", length(rp), "| of these in the 51:", length(intersect(rp, p51)),
"| the 16 modelled that are RPPA patients:", length(intersect(p16, rp)), "\n")
}
}
pm <- fread("01_RawData/Breast/brca_prot_meta.csv")
cat("\n76-patient matrix, Group column:\n"); print(table(pm$Group, useNA = "ifany"))
})
sec("5. PRAD RNA INPUT: which scale was saved, and is a log-scale version available?", {
cat("files in 02_ProcessedData mentioning vst, log or rna:\n")
print(grep("vst|log|rna", list.files("02_ProcessedData", recursive = TRUE), ignore.case = TRUE, value = TRUE))
f <- "02_ProcessedData/Prostate/prad_rna_normalized.csv"
if (ex(f)) {
x <- fread(f, nrows = 3000)
num <- x[, sapply(x, is.numeric), with = FALSE]
cat("first 3000 genes of", f, ": range", range(unlist(num), na.rm = TRUE), "| median", median(unlist(num), na.rm = TRUE), "\n")
}
e <- new.env(); load(file.path(b, "PRAD/prep_objects.RData"), envir = e)
r <- as.matrix(e$A$X$rna)
cat("PRAD DIABLO/MOFA RNA block: quantiles of all values\n"); print(quantile(r, c(0, .25, .5, .75, .9, .99, .999, 1)))
v <- apply(r, 1, var); top <- order(v, decreasing = TRUE)[1:5000]
cat("share of total variance held by the 20 highest-variance genes:", round(sum(sort(v, decreasing = TRUE)[1:20]) / sum(v), 3), "\n")
cat("the 10 highest-variance genes:", paste(head(rownames(r)[order(v, decreasing = TRUE)], 10), collapse = ", "), "\n")
})
sec("6. CLASSICAL PROSTATE MARKERS WITH WHOLE-SYMBOL MATCHING (writes a new file)", {
suppressPackageStartupMessages({ library(minfi); library(IlluminaHumanMethylation450kanno.ilmn12.hg19) })
dmp <- fread("03_Results/Epigenomics/prad_DMPs.csv")
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
gn  <- ann[match(dmp$probe_id, rownames(ann)), "UCSC_RefGene_Name"]
sym <- strsplit(as.character(gn), ";")
out <- rbindlist(lapply(c("GSTP1", "RASSF1", "APC", "PTGS2"), function(g) {
h <- vapply(sym, function(x) g %in% x, logical(1))
d <- dmp[h, .(probe_id, delta_beta, logFC_Mvalue = logFC, adj.P.Val)]
d[, `:=`(gene = gn[h], marker = g)]
d[order(adj.P.Val)]
}))
fwrite(out, "classic_marker_validation_exact.csv")
print(out[, .(CpGs = .N, smallest_adjP = min(adj.P.Val), higher_in_cancer = sum(delta_beta > 0), lower_in_cancer = sum(delta_beta < 0),
min_delta = min(delta_beta), max_delta = max(delta_beta)), by = marker])
cat("written: classic_marker_validation_exact.csv\n")
})
cat("\n\nDONE. Upload needed_outputs_2.txt.\n")
while (sink.number() > 0) sink()
cat("DONE. Output written to", file.path(getwd(), "needed_outputs_2.txt"), "\n")
