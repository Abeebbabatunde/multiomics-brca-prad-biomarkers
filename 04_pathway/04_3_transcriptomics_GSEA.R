# ============================================================
# STEP 4.3 -- Transcriptomics GSEA (KEGG, GO-BP), BRCA and PRAD
# Inputs : 03_Results/Transcriptomics/brca_DEGs_all.csv (first column `...1`, Ensembl IDs with version; `stat`)
#          03_Results/Transcriptomics/prad_DEGs_full_results.csv (no `stat`; ranked by sign(log2FC) * -log10(p))
# Outputs: 03_Results/Transcriptomics/{BRCA,PRAD}_transcriptomics_{KEGG,GO_BP}_GSEA.rds
# Console: BRCA 44,607 genes -> 29,103 mapped rows -> 28,985 ranked Entrez IDs
#          PRAD 43,076 genes -> 28,393 mapped rows -> 28,269 ranked Entrez IDs; ties in 1.43% (PRAD) / 1.42% (BRCA)
#          Result rows (Step 5.0 inventory): BRCA GO 2082, BRCA KEGG 167, PRAD GO 1405, PRAD KEGG 142
# No set.seed in this version (the earlier 04_0 first-run script had set.seed(1234)). GSEA uses permutations: rerun values may differ slightly.
# ============================================================
library(readr); library(dplyr); library(clusterProfiler); library(org.Hs.eg.db)

# ---- BRCA ----
brca_deg <- read_csv("03_Results/Transcriptomics/brca_DEGs_all.csv", show_col_types = FALSE)
brca_deg <- brca_deg %>%
  mutate(ensembl_clean = sub("\\..*$", "", `...1`)) %>%
  filter(!is.na(stat), !is.na(ensembl_clean)) %>%
  distinct(ensembl_clean, .keep_all = TRUE)
brca_map <- bitr(brca_deg$ensembl_clean, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
brca_deg <- inner_join(brca_deg, brca_map, by = c("ensembl_clean" = "ENSEMBL"))
brca_ranked <- brca_deg %>% arrange(desc(stat)) %>% distinct(ENTREZID, .keep_all = TRUE)
brca_gene_list <- setNames(brca_ranked$stat, brca_ranked$ENTREZID)

brca_deg_kegg <- gseKEGG(geneList = brca_gene_list, organism = "hsa", minGSSize = 10, maxGSSize = 500,
                         pvalueCutoff = 0.05, eps = 0)
saveRDS(brca_deg_kegg, "03_Results/Transcriptomics/BRCA_transcriptomics_KEGG_GSEA.rds")
brca_deg_go <- gseGO(geneList = brca_gene_list, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                     minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, eps = 0)
saveRDS(brca_deg_go, "03_Results/Transcriptomics/BRCA_transcriptomics_GO_BP_GSEA.rds")

# ---- PRAD ----
prad_deg <- read_csv("03_Results/Transcriptomics/prad_DEGs_full_results.csv", show_col_types = FALSE)
prad_deg <- prad_deg %>%
  mutate(ensembl_clean = sub("\\..*$", "", ensembl_gene_id),
         rank_metric = sign(log2FoldChange) * -log10(pmax(pvalue, 1e-300))) %>%
  filter(!is.na(rank_metric), is.finite(rank_metric)) %>%
  distinct(ensembl_clean, .keep_all = TRUE)
prad_map <- bitr(prad_deg$ensembl_clean, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
prad_deg <- inner_join(prad_deg, prad_map, by = c("ensembl_clean" = "ENSEMBL"))
prad_ranked <- prad_deg %>% arrange(desc(rank_metric)) %>% distinct(ENTREZID, .keep_all = TRUE)
prad_gene_list <- setNames(prad_ranked$rank_metric, prad_ranked$ENTREZID)

prad_deg_kegg <- gseKEGG(geneList = prad_gene_list, organism = "hsa", minGSSize = 10, maxGSSize = 500,
                         pvalueCutoff = 0.05, eps = 0)
saveRDS(prad_deg_kegg, "03_Results/Transcriptomics/PRAD_transcriptomics_KEGG_GSEA.rds")
prad_deg_go <- gseGO(geneList = prad_gene_list, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                     minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, eps = 0)
saveRDS(prad_deg_go, "03_Results/Transcriptomics/PRAD_transcriptomics_GO_BP_GSEA.rds")
