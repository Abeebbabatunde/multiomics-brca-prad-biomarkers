# ============================================================
# STEP 4.4 -- BRCA proteomics (RPPA) ORA, KEGG + GO-BP, improved (direct SYMBOL + one-to-one ALIAS) mapping
# Input : 03_Results/Proteomics/brca_DEPs.csv ; antibody names stripped of phospho suffix with sub("_p[STY][0-9]+.*$","")
# Output: BRCA_proteomics_KEGG_ORA_corrected.rds, BRCA_proteomics_GO_BP_ORA_corrected.rds (03_Results/Proteomics/)
# Console: 402 cleaned targets; direct SYMBOL 162; unmapped 240; ALIAS-recoverable 90 (124 records),
#          one-to-one 72 / ambiguous 18 (e.g. PAR, PR, RIP, DRP1); mapped 234; unmapped 168;
#          background 234, significant 69 (29.5%); KEGG 0 terms, GO-BP 0 terms.
# Superseded first attempt (direct SYMBOL only): 162 background, 49 significant -> BRCA_proteomics_KEGG_ORA.rds (0 terms).
# Note: the GO-BP call in the "corrected" console run had no qvalueCutoff argument; the 04_0-style script set qvalueCutoff = 0.05.
# KEGG: online lookup (rest.kegg.jp); retrieval date not recorded.
# ============================================================
library(readr); library(dplyr); library(clusterProfiler); library(org.Hs.eg.db)

brca_dep_check <- read_csv("03_Results/Proteomics/brca_DEPs.csv", show_col_types = FALSE) %>%
  mutate(gene_symbol = sub("_p[STY][0-9]+.*$", "", peptide_target)) %>%
  arrange(p_value) %>%
  distinct(gene_symbol, .keep_all = TRUE)

brca_dep_check_map <- suppressWarnings(bitr(brca_dep_check$gene_symbol, fromType = "SYMBOL",
                                            toType = "ENTREZID", OrgDb = org.Hs.eg.db))
brca_unmapped <- setdiff(brca_dep_check$gene_symbol, brca_dep_check_map$SYMBOL)

brca_alias_map <- suppressWarnings(bitr(brca_unmapped, fromType = "ALIAS", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
brca_alias_one_to_one <- brca_alias_map %>%
  group_by(ALIAS) %>% filter(n_distinct(ENTREZID) == 1) %>%
  summarise(ENTREZID = unique(ENTREZID), .groups = "drop")

brca_dep_improved <- brca_dep_check %>%
  left_join(brca_dep_check_map %>% distinct(SYMBOL, ENTREZID), by = c("gene_symbol" = "SYMBOL")) %>%
  rename(ENTREZID_direct = ENTREZID) %>%
  left_join(brca_alias_one_to_one, by = c("gene_symbol" = "ALIAS")) %>%
  rename(ENTREZID_alias = ENTREZID) %>%
  mutate(ENTREZID_final = ifelse(!is.na(ENTREZID_direct), ENTREZID_direct, ENTREZID_alias))

brca_dep_mapped <- brca_dep_improved %>% filter(!is.na(ENTREZID_final))
brca_dep_universe_improved <- unique(brca_dep_mapped$ENTREZID_final)                       # 234
brca_dep_sig_improved <- unique(brca_dep_mapped$ENTREZID_final[brca_dep_mapped$significant == TRUE])  # 69

brca_dep_kegg_corrected <- enrichKEGG(gene = brca_dep_sig_improved, universe = brca_dep_universe_improved,
                                      organism = "hsa", minGSSize = 5, maxGSSize = 500, pvalueCutoff = 0.05)
saveRDS(brca_dep_kegg_corrected, "03_Results/Proteomics/BRCA_proteomics_KEGG_ORA_corrected.rds")

brca_dep_go_corrected <- enrichGO(gene = brca_dep_sig_improved, universe = brca_dep_universe_improved,
                                  OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                                  minGSSize = 5, maxGSSize = 500, pvalueCutoff = 0.05)
saveRDS(brca_dep_go_corrected, "03_Results/Proteomics/BRCA_proteomics_GO_BP_ORA_corrected.rds")
