## 05_5  PRAD classical methylation markers -- script AS RUN (produced classic_marker_validation.csv, 89 rows)
## WARNING: grepl(g, dmp$gene) is an unanchored substring match: "RASSF1" also matches RASSF10 and "APC" matches
## ANAPC2, APCDD1, SNAPC2 etc. Whole-symbol versions:
## 05_5b (original table, classic_marker_validation_exact.csv) and 05_5c (corrected 79,679-probe table,
## classic_marker_validation_exact_corrected.csv), both recovered from RStudio history.
setwd("C:/Users/USER/Documents/MyThesis_MultiOmics")
suppressPackageStartupMessages({ library(data.table); library(minfi); library(IlluminaHumanMethylation450kanno.ilmn12.hg19) })
dmp <- fread("03_Results/Epigenomics/prad_DMPs.csv")
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
dmp$gene <- ann[match(dmp$probe_id, rownames(ann)), "UCSC_RefGene_Name"]
classic <- c("GSTP1","RASSF1","APC","PTGS2")
for (g in classic) {
  hits <- dmp[grepl(g, dmp$gene), c("probe_id","logFC","adj.P.Val","gene")]
  cat("\n", g, ":\n"); print(hits)
}
