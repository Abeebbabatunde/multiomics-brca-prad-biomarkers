# ============================================================
# STEP 4.5 -- PRAD proteomics (mass spectrometry) ORA, nominal P < 0.05, one-to-one SYMBOL->Entrez
# Input : 03_Results/Proteomics/prad_DEPs_all.tsv (5,412 proteins)
# Output: PRAD_proteomics_{clean_mapped,background_Entrez,FDR_sig_Entrez,nominal_sig_Entrez}.rds,
#         PRAD_proteomics_KEGG_ORA_nominal.rds, PRAD_proteomics_GO_BP_ORA_nominal.rds
# Console: 5,412 rows; 5,338 mapped records / 5,337 symbols; HBD maps to 3045 and 100187828 -> removed;
#          clean one-to-one 5,336; FDR < 0.05: 4; nominal P < 0.05: 718 (13.5%)
#          GO-BP: 1 term (cell-substrate adhesion, 39/696, P 5.68e-06, adj P 0.0285)
#          KEGG : 1 term (hsa04820 Cytoskeleton in muscle cells, 31/420, P 9.78e-05, adj P 0.0329)
# The nominal threshold is exploratory; FDR 4 genes is too few for ORA.
# ============================================================
library(readr); library(dplyr); library(clusterProfiler); library(org.Hs.eg.db)

prad_dep <- read_delim("03_Results/Proteomics/prad_DEPs_all.tsv", delim = "\t", show_col_types = FALSE) %>%
  filter(!is.na(GeneName), GeneName != "") %>%
  distinct(GeneName, .keep_all = TRUE)

prad_dep_map <- suppressWarnings(bitr(prad_dep$GeneName, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
prad_dep_map_clean <- prad_dep_map %>%
  group_by(SYMBOL) %>% filter(n_distinct(ENTREZID) == 1) %>%
  summarise(ENTREZID = unique(ENTREZID), .groups = "drop")                     # 5336

prad_dep_clean <- prad_dep %>% inner_join(prad_dep_map_clean, by = c("GeneName" = "SYMBOL"))

prad_dep_universe_clean    <- unique(prad_dep_clean$ENTREZID)
prad_dep_sig_fdr_clean     <- unique(prad_dep_clean$ENTREZID[!is.na(prad_dep_clean$FDR) & prad_dep_clean$FDR < 0.05])
prad_dep_sig_nominal_clean <- unique(prad_dep_clean$ENTREZID[!is.na(prad_dep_clean$p_value) & prad_dep_clean$p_value < 0.05])

saveRDS(prad_dep_clean,               "03_Results/Proteomics/PRAD_proteomics_clean_mapped.rds")
saveRDS(prad_dep_universe_clean,      "03_Results/Proteomics/PRAD_proteomics_background_Entrez.rds")
saveRDS(prad_dep_sig_fdr_clean,       "03_Results/Proteomics/PRAD_proteomics_FDR_sig_Entrez.rds")
saveRDS(prad_dep_sig_nominal_clean,   "03_Results/Proteomics/PRAD_proteomics_nominal_sig_Entrez.rds")

prad_prot_kegg_nominal <- enrichKEGG(gene = prad_dep_sig_nominal_clean, universe = prad_dep_universe_clean,
                                     organism = "hsa", minGSSize = 5, maxGSSize = 500, pvalueCutoff = 0.05)
saveRDS(prad_prot_kegg_nominal, "03_Results/Proteomics/PRAD_proteomics_KEGG_ORA_nominal.rds")

prad_prot_go_nominal <- enrichGO(gene = prad_dep_sig_nominal_clean, universe = prad_dep_universe_clean,
                                 OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                                 minGSSize = 5, maxGSSize = 500, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
saveRDS(prad_prot_go_nominal, "03_Results/Proteomics/PRAD_proteomics_GO_BP_ORA_nominal.rds")
