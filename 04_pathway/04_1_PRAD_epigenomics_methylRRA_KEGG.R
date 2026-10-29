# ============================================================
# STEP 4.1 -- PRAD epigenomics: methylRRA (ORA, KEGG)
# Input : 03_Results/Epigenomics/prad_DMPs_full_results.csv (485,577 probes; loaded as `prad_meth`)
# Output: 03_Results/Epigenomics/prad_methylRRA_KEGG.rds
# Console: 64,426 probes have no P value -> dropped; 421,151 submitted;
#          "213 gene sets are being tested", "19066 genes are significant under cut-off 0.05"
# NOTE: the first attempt (no NA filter) failed: "missing value where TRUE/FALSE needed".
# NOT SUPPLIED: the BRCA methylRRA (KEGG/GO) and PRAD methylRRA GO writers; no such RDS is in the Step 5.0 inventory.
# Prerequisites: library(methylGSA); library(readr); ARRAY_TYPE <- "450K"
# ============================================================
library(methylGSA); library(readr)
ARRAY_TYPE <- "450K"
prad_meth <- read_csv("03_Results/Epigenomics/prad_DMPs_full_results.csv", show_col_types = FALSE)

prad_meth_clean <- prad_meth[
  !is.na(prad_meth$probe_id) & !is.na(prad_meth$P.Value) & is.finite(prad_meth$P.Value), ]
prad_meth_pvals <- setNames(prad_meth_clean$P.Value, prad_meth_clean$probe_id)

# checks (all printed 0 / expected values): length 421151; NA 0; ==0 0; <0 0; >1 0
prad_dmp_kegg <- methylRRA(cpg.pval = prad_meth_pvals, method = "ORA", array.type = ARRAY_TYPE,
                           GS.type = "KEGG", minsize = 10, maxsize = 500)
saveRDS(prad_dmp_kegg, "03_Results/Epigenomics/prad_methylRRA_KEGG.rds")
