## ARCHIVE -- NOT CONFIRMED RUN AS WRITTEN. No console log accompanies this file.
## Part 1 is a BRCA methylation limma template: it needs `brca_clinical$group` (Tumor/Normal), which the script itself says to
##   replace, and reads 01_RawData/brca_methylation_raw.rds. The script that produced the reported results is
##   ../03_2_BRCA_methylation_limma_AS_RUN.R (305,553 probes tested, 18,565 DMPs). This template calls its DMP cut-off |logFC| >= 0.10.
## Part 2 is an early draft of the Step 4 methylRRA (PRAD KEGG) and BRCA transcriptomics GSEA code; it assumes `prad_meth` and
##   ARRAY_TYPE already exist. The scripts used are in ../../04_pathway/. The closing sentence of the paste ("Important: this
##   combines the two analyses...") is a note that came with the code and is kept as the last comment line.
# ============================================================
# BRCA METHYLATION ANALYSIS
# Raw Data → QC/Filtering → Differential Methylation
# → Significant DMPs → Save Results
# ============================================================


# ============================================================
# 1. LOAD REQUIRED PACKAGES
# ============================================================

library(limma)
library(dplyr)
library(tibble)


# ============================================================
# 2. LOCATE BRCA METHYLATION FILES
# ============================================================

list.files(
  "01_RawData",
  pattern = "brca.*methyl|methyl.*brca",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)


# ============================================================
# 3. LOAD RAW BRCA METHYLATION OBJECT
# ============================================================

brca_methylation_raw <- readRDS(
  "01_RawData/brca_methylation_raw.rds"
)

dim(brca_methylation_raw)


# ============================================================
# 4. CONVERT TO MATRIX
# ============================================================

brca_methylation <- as.matrix(
  brca_methylation_raw
)

dim(brca_methylation)


# ============================================================
# 5. RECORD INITIAL NUMBER OF PROBES
# ============================================================

brca_n_start <- nrow(
  brca_methylation
)

cat(
  "Starting number of probes:",
  brca_n_start,
  "\n"
)


# ============================================================
# 6. FILTER PROBES WITH EXCESSIVE MISSING VALUES
#    Criterion: ≤20% missing values per probe
# ============================================================

brca_na_prop <- rowMeans(
  is.na(brca_methylation)
)

brca_methylation_filt <- brca_methylation[
  brca_na_prop <= 0.20,
  ,
  drop = FALSE
]

brca_n_after_na <- nrow(
  brca_methylation_filt
)

cat(
  "After NA filtering:",
  brca_n_after_na,
  "\n"
)

cat(
  "Probes removed by NA filtering:",
  brca_n_start - brca_n_after_na,
  "\n"
)


# ============================================================
# 7. REMOVE PROBES WITH REMAINING MISSING VALUES
# ============================================================

brca_methylation_filt <- brca_methylation_filt[
  complete.cases(brca_methylation_filt),
  ,
  drop = FALSE
]

brca_n_after_complete <- nrow(
  brca_methylation_filt
)

cat(
  "After complete-case filtering:",
  brca_n_after_complete,
  "\n"
)


# ============================================================
# 8. REMOVE ZERO-VARIANCE PROBES
# ============================================================

brca_probe_sd <- apply(
  brca_methylation_filt,
  1,
  sd,
  na.rm = TRUE
)

brca_methylation_filt <- brca_methylation_filt[
  brca_probe_sd > 0,
  ,
  drop = FALSE
]

brca_n_after_variance <- nrow(
  brca_methylation_filt
)

cat(
  "After variance filtering:",
  brca_n_after_variance,
  "\n"
)


# ============================================================
# 9. LOAD / PREPARE BRCA CLINICAL METADATA
# ============================================================

# Replace this with the actual clinical metadata object
# if it already exists in your workspace.

# Example:
#
# brca_clinical <- readRDS(
#   "01_RawData/brca_clinical.rds"
# )


# ============================================================
# 10. DEFINE BRCA GROUP
# ============================================================

# IMPORTANT:
# Replace "group" with the actual column containing
# Tumor/Normal information.

brca_group <- factor(
  brca_clinical$group,
  levels = c(
    "Normal",
    "Tumor"
  )
)

table(brca_group)


# ============================================================
# 11. CHECK SAMPLE NUMBER AGREEMENT
# ============================================================

cat(
  "Number of methylation samples:",
  ncol(brca_methylation_filt),
  "\n"
)

cat(
  "Number of clinical samples:",
  length(brca_group),
  "\n"
)


# ============================================================
# 12. CREATE DESIGN MATRIX
# ============================================================

brca_design <- model.matrix(
  ~ brca_group
)

colnames(brca_design)


# ============================================================
# 13. FIT LIMMA MODEL
# ============================================================

brca_fit <- lmFit(
  brca_methylation_filt,
  brca_design
)


# ============================================================
# 14. EMPIRICAL BAYES MODERATION
# ============================================================

brca_fit <- eBayes(
  brca_fit
)


# ============================================================
# 15. EXTRACT COMPLETE DIFFERENTIAL
#    METHYLATION RESULTS
# ============================================================

brca_diffmeth_results <- topTable(
  brca_fit,
  coef = "brca_groupTumor",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)


# ============================================================
# 16. ADD PROBE IDs
# ============================================================

brca_diffmeth_results <- brca_diffmeth_results %>%
  rownames_to_column(
    "Probe_ID"
  )


# ============================================================
# 17. DEFINE SIGNIFICANT DIFFERENTIALLY
#     METHYLATED PROBES
#
#     Criteria:
#     Adjusted P < 0.05
#     Absolute logFC ≥ 0.10
# ============================================================

brca_dmp_sig <- brca_diffmeth_results %>%
  filter(
    adj.P.Val < 0.05,
    abs(logFC) >= 0.10
  )


# ============================================================
# 18. REPORT RESULTS
# ============================================================

cat(
  "\n====================================\n"
)

cat(
  "BRCA METHYLATION ANALYSIS SUMMARY\n"
)

cat(
  "====================================\n"
)

cat(
  "Starting probes:",
  brca_n_start,
  "\n"
)

cat(
  "After NA filtering:",
  brca_n_after_na,
  "\n"
)

cat(
  "After complete-case filtering:",
  brca_n_after_complete,
  "\n"
)

cat(
  "After variance filtering:",
  brca_n_after_variance,
  "\n"
)

cat(
  "Total probes tested:",
  nrow(brca_diffmeth_results),
  "\n"
)

cat(
  "Significant DMPs:",
  nrow(brca_dmp_sig),
  "\n"
)


# ============================================================
# 19. SAVE COMPLETE DIFFERENTIAL METHYLATION RESULTS
# ============================================================

saveRDS(
  brca_diffmeth_results,
  "02_ProcessedData/breast_methylation_diffmeth_results.rds"
)


# ============================================================
# 20. SAVE SIGNIFICANT DMPs
# ============================================================

saveRDS(
  brca_dmp_sig,
  "02_ProcessedData/breast_methylation_significant_DMPs.rds"
)


# ============================================================
# 21. CREATE EXPLICIT TESTED-PROBE UNIVERSE
#     FOR DOWNSTREAM ANALYSIS / methylRRA
# ============================================================

brca_methylation_universe <- rownames(
  brca_methylation_filt
)

saveRDS(
  brca_methylation_universe,
  "02_ProcessedData/breast_methylation_tested_probe_universe.rds"
)


# ============================================================
# 22. CHECK FINAL FILES
# ============================================================

file.exists(
  "02_ProcessedData/breast_methylation_diffmeth_results.rds"
)

file.exists(
  "02_ProcessedData/breast_methylation_significant_DMPs.rds"
)

file.exists(
  "02_ProcessedData/breast_methylation_tested_probe_universe.rds"
)


# ============================================================
# 23. FINAL DIMENSION CHECK
# ============================================================

cat(
  "\nFinal differential methylation results:",
  nrow(brca_diffmeth_results),
  "probes\n"
)

cat(
  "Final significant DMPs:",
  nrow(brca_dmp_sig),
  "\n"
)

cat(
  "Final tested-probe universe:",
  length(brca_methylation_universe),
  "probes\n"
)








## ============================================================
## MULTI-OMICS PATHWAY ANALYSIS
## PRAD EPIGENOMICS + BRCA TRANSCRIPTOMICS
## ============================================================


## ============================================================
## 1. PRAD EPIGENOMICS — methylRRA KEGG
## ============================================================

prad_meth_clean <- prad_meth[
  !is.na(prad_meth$probe_id) &
    !is.na(prad_meth$P.Value) &
    is.finite(prad_meth$P.Value),
]

prad_meth_pvals <- setNames(
  prad_meth_clean$P.Value,
  prad_meth_clean$probe_id
)

prad_dmp_kegg <- methylRRA(
  cpg.pval = prad_meth_pvals,
  method = "ORA",
  array.type = ARRAY_TYPE,
  GS.type = "KEGG",
  minsize = 10,
  maxsize = 500
)

saveRDS(
  prad_dmp_kegg,
  "03_Results/Epigenomics/prad_methylRRA_KEGG.rds"
)


## ============================================================
## 2. BRCA TRANSCRIPTOMICS — Prepare ranked gene list
## ============================================================

library(readr)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

brca_deg <- read_csv(
  "03_Results/Transcriptomics/brca_DEGs_all.csv",
  show_col_types = FALSE
)

brca_deg <- brca_deg %>%
  mutate(
    ensembl_clean = sub("\\..*$", "", `...1`)
  ) %>%
  filter(
    !is.na(stat),
    !is.na(ensembl_clean)
  ) %>%
  distinct(
    ensembl_clean,
    .keep_all = TRUE
  )

brca_map <- bitr(
  brca_deg$ensembl_clean,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

brca_deg <- inner_join(
  brca_deg,
  brca_map,
  by = c("ensembl_clean" = "ENSEMBL")
)

brca_ranked <- brca_deg %>%
  arrange(desc(stat)) %>%
  distinct(
    ENTREZID,
    .keep_all = TRUE
  )

brca_gene_list <- setNames(
  brca_ranked$stat,
  brca_ranked$ENTREZID
)


## ============================================================
## 3. BRCA TRANSCRIPTOMICS — KEGG GSEA
## ============================================================

brca_deg_kegg <- gseKEGG(
  geneList = brca_gene_list,
  organism = "hsa",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  eps = 0
)

saveRDS(
  brca_deg_kegg,
  "03_Results/Transcriptomics/BRCA_transcriptomics_KEGG_GSEA.rds"
)


## ============================================================
## 4. BRCA TRANSCRIPTOMICS — GO Biological Process GSEA
## ============================================================

brca_deg_go <- gseGO(
  geneList = brca_gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  eps = 0
)

saveRDS(
  brca_deg_go,
  "03_Results/Transcriptomics/BRCA_transcriptomics_GO_BP_GSEA.rds"
)


## ============================================================
## END
## ============================================================

## Important: this combines the two analyses we have actually established. The PRAD section assumes "prad_meth", "ARRAY_TYPE", and the "methylRRA" package/function are already loaded, just as in your existing workflow.