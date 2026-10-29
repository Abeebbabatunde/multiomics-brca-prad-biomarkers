## =============================================================================
## STEP 4.0 — Individual Omics Pathway Analysis (FINAL, wired to real files)
## Working directory assumed: MyThesis_MultiOmics/ (run setwd() first if needed)
## =============================================================================
##
## Source files confirmed from your project:
##   BRCA epigenomics : results/step2_2_methylation_diff/breast_methylation_diffmeth_results.rds
##                       (305,553 probes, FULL — post-QC-filtered 450K array; NOT the
##                        18,565-row 03_Results/Epigenomics/brca_DMPs.csv, which is
##                        pre-filtered to significant hits only and is NOT suitable
##                        input for methylRRA)
##   PRAD epigenomics : 03_Results/Epigenomics/prad_DMPs_full_results.csv (485,577 probes, FULL)
##   BRCA transcriptomics : 03_Results/Transcriptomics/brca_DEGs_all.csv (44,607 genes, FULL)
##   PRAD transcriptomics : 03_Results/Transcriptomics/prad_DEGs_full_results.csv (43,076 genes, FULL)
##   BRCA proteomics  : 03_Results/Proteomics/brca_DEPs.csv (457 antibodies, RPPA panel)
##   PRAD proteomics  : 03_Results/Proteomics/prad_DEPs_all.tsv (5,412 proteins, FULL, mass-spec)
##
## FLAGGED ASSUMPTIONS — confirm before trusting results:
##   1. Array type set to "450K" for both BRCA and PRAD methylation (485,577 = exact
##      450K probe count; BRCA's 305,553 is a QC-filtered subset of the same array).
##      Change ARRAY_TYPE below if either was actually run on EPIC/850K.
##   2. PRAD DEGs have no "stat" column (BRCA's does) — GSEA rank built instead from
##      sign(log2FoldChange) * -log10(pvalue). Less ideal than a Wald/LRT stat but
##      the closest available substitute.
##   3. BRCA proteomics is an RPPA panel (457 antibodies, NOT full proteome) — pathway
##      enrichment here has inherently low power/resolution vs PRAD's mass-spec data
##      (5,412 proteins). Antibody names are cleaned of phospho-site suffixes
##      (e.g. "CDK1_pT14" -> "CDK1") before symbol mapping; multiple antibodies per
##      gene are deduplicated by keeping the most significant.
##   4. PRAD proteomics: only 4 proteins pass FDR < 0.05 (prad_DEPs_significant.tsv) —
##      too few for meaningful ORA. Script below defaults to NOMINAL p < 0.05 for the
##      proteomics significance cut specifically (documented at that step) to get a
##      workable gene list; switch back to FDR if you disagree.
## =============================================================================

required_pkgs <- c("methylGSA", "missMethyl", "clusterProfiler",
                    "org.Hs.eg.db", "fgsea", "ReactomePA", "dplyr", "readr")
# BiocManager::install(required_pkgs)  # run once if any are missing

library(methylGSA)
library(missMethyl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(dplyr)
library(readr)

set.seed(1234)
ARRAY_TYPE <- "450K"   # change to "EPIC" if either methylation dataset used that platform

out_dir <- "04_Pathway_Analysis"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (sub in c("Epigenomics", "Transcriptomics", "Proteomics")) {
  dir.create(file.path(out_dir, sub), showWarnings = FALSE)
}

save_enrich <- function(res_obj, sub, label, is_gsea = FALSE) {
  if (is.null(res_obj)) { cat("Skipped (NULL result):", label, "\n"); return(invisible(NULL)) }
  df <- if (inherits(res_obj, "data.frame")) res_obj else as.data.frame(res_obj)
  write_csv(df, file.path(out_dir, sub, paste0(label, "_results.csv")))
  cat("Saved:", label, "-", nrow(df), "pathways\n")
  invisible(df)
}

## =============================================================================
## 1. BRCA EPIGENOMICS — full 450K array (305,553 probes) via methylRRA
## =============================================================================
brca_meth <- readRDS("results/step2_2_methylation_diff/breast_methylation_diffmeth_results.rds")
brca_meth_pvals <- setNames(brca_meth$P.Value, brca_meth$Probe)

brca_dmp_kegg <- methylRRA(
  cpg.pval = brca_meth_pvals, method = "ORA",
  array.type = ARRAY_TYPE, GS.type = "KEGG", minsize = 10, maxsize = 500
)
save_enrich(brca_dmp_kegg, "Epigenomics", "BRCA_epigenomics_KEGG")

brca_dmp_go <- methylRRA(
  cpg.pval = brca_meth_pvals, method = "ORA",
  array.type = ARRAY_TYPE, GS.type = "GO", minsize = 10, maxsize = 500
)
save_enrich(brca_dmp_go, "Epigenomics", "BRCA_epigenomics_GO")

## =============================================================================
## 2. PRAD EPIGENOMICS — full 450K array (485,577 probes) via methylRRA
## =============================================================================
prad_meth <- read_csv("03_Results/Epigenomics/prad_DMPs_full_results.csv", show_col_types = FALSE)
prad_meth_pvals <- setNames(prad_meth$P.Value, prad_meth$probe_id)

prad_dmp_kegg <- methylRRA(
  cpg.pval = prad_meth_pvals, method = "ORA",
  array.type = ARRAY_TYPE, GS.type = "KEGG", minsize = 10, maxsize = 500
)
save_enrich(prad_dmp_kegg, "Epigenomics", "PRAD_epigenomics_KEGG")

prad_dmp_go <- methylRRA(
  cpg.pval = prad_meth_pvals, method = "ORA",
  array.type = ARRAY_TYPE, GS.type = "GO", minsize = 10, maxsize = 500
)
save_enrich(prad_dmp_go, "Epigenomics", "PRAD_epigenomics_GO")

## =============================================================================
## 3. BRCA TRANSCRIPTOMICS — full DESeq2 results (44,607 genes) via GSEA
## =============================================================================
brca_deg <- read_csv("03_Results/Transcriptomics/brca_DEGs_all.csv", show_col_types = FALSE)

brca_deg <- brca_deg %>%
  mutate(ensembl_clean = sub("\\..*$", "", X)) %>%     # strip Ensembl version suffix
  filter(!is.na(stat), !is.na(ensembl_clean)) %>%
  distinct(ensembl_clean, .keep_all = TRUE)

brca_map <- bitr(brca_deg$ensembl_clean, "ENSEMBL", "ENTREZID", org.Hs.eg.db)
brca_deg <- inner_join(brca_deg, brca_map, by = c("ensembl_clean" = "ENSEMBL"))

brca_ranked <- brca_deg %>% arrange(desc(stat)) %>% distinct(ENTREZID, .keep_all = TRUE)
brca_gene_list <- setNames(brca_ranked$stat, brca_ranked$ENTREZID)

brca_deg_kegg <- gseKEGG(geneList = brca_gene_list, organism = "hsa",
                          minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, eps = 0)
save_enrich(brca_deg_kegg, "Transcriptomics", "BRCA_transcriptomics_KEGG_GSEA", is_gsea = TRUE)

brca_deg_go <- gseGO(geneList = brca_gene_list, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                      ont = "BP", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, eps = 0)
save_enrich(brca_deg_go, "Transcriptomics", "BRCA_transcriptomics_GO_BP_GSEA", is_gsea = TRUE)

## =============================================================================
## 4. PRAD TRANSCRIPTOMICS — full DESeq2 results (43,076 genes) via GSEA
##    NOTE: no "stat" column available — rank built from sign(log2FC) * -log10(pvalue)
## =============================================================================
prad_deg <- read_csv("03_Results/Transcriptomics/prad_DEGs_full_results.csv", show_col_types = FALSE)

prad_deg <- prad_deg %>%
  mutate(ensembl_clean = sub("\\..*$", "", ensembl_gene_id),
         rank_metric = sign(log2FoldChange) * -log10(pmax(pvalue, 1e-300))) %>%
  filter(!is.na(rank_metric), is.finite(rank_metric)) %>%
  distinct(ensembl_clean, .keep_all = TRUE)

prad_map <- bitr(prad_deg$ensembl_clean, "ENSEMBL", "ENTREZID", org.Hs.eg.db)
prad_deg <- inner_join(prad_deg, prad_map, by = c("ensembl_clean" = "ENSEMBL"))

prad_ranked <- prad_deg %>% arrange(desc(rank_metric)) %>% distinct(ENTREZID, .keep_all = TRUE)
prad_gene_list <- setNames(prad_ranked$rank_metric, prad_ranked$ENTREZID)

prad_deg_kegg <- gseKEGG(geneList = prad_gene_list, organism = "hsa",
                          minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, eps = 0)
save_enrich(prad_deg_kegg, "Transcriptomics", "PRAD_transcriptomics_KEGG_GSEA", is_gsea = TRUE)

prad_deg_go <- gseGO(geneList = prad_gene_list, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                      ont = "BP", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, eps = 0)
save_enrich(prad_deg_go, "Transcriptomics", "PRAD_transcriptomics_GO_BP_GSEA", is_gsea = TRUE)

## =============================================================================
## 5. BRCA PROTEOMICS — RPPA panel (457 antibodies) via ORA, custom background
## =============================================================================
brca_dep <- read_csv("03_Results/Proteomics/brca_DEPs.csv", show_col_types = FALSE)

# Clean antibody names: strip phospho-site suffixes (e.g. "CDK1_pT14" -> "CDK1")
brca_dep <- brca_dep %>%
  mutate(gene_symbol = sub("_p[STY][0-9]+.*$", "", peptide_target)) %>%
  arrange(p_value) %>%
  distinct(gene_symbol, .keep_all = TRUE)   # keep most significant antibody per gene

brca_dep_map <- suppressWarnings(bitr(brca_dep$gene_symbol, "SYMBOL", "ENTREZID", org.Hs.eg.db))
cat("BRCA RPPA: mapped", nrow(brca_dep_map), "of", nrow(brca_dep), "antibody targets to Entrez IDs\n")

brca_dep <- inner_join(brca_dep, brca_dep_map, by = c("gene_symbol" = "SYMBOL"))
brca_dep_universe <- unique(brca_dep$ENTREZID)
brca_dep_sig <- unique(brca_dep$ENTREZID[brca_dep$significant == TRUE])

cat("BRCA RPPA: ", length(brca_dep_sig), "significant of", length(brca_dep_universe),
    "mapped background (RPPA is a targeted 457-antibody panel — treat results as exploratory)\n")

brca_dep_kegg <- if (length(brca_dep_sig) >= 3) {
  enrichKEGG(gene = brca_dep_sig, universe = brca_dep_universe, organism = "hsa",
             minGSSize = 5, maxGSSize = 500, pvalueCutoff = 0.05)
} else NULL
save_enrich(brca_dep_kegg, "Proteomics", "BRCA_proteomics_KEGG_ORA")

brca_dep_go <- if (length(brca_dep_sig) >= 3) {
  enrichGO(gene = brca_dep_sig, universe = brca_dep_universe, OrgDb = org.Hs.eg.db,
            keyType = "ENTREZID", ont = "BP", minGSSize = 5, maxGSSize = 500, pvalueCutoff = 0.05)
} else NULL
save_enrich(brca_dep_go, "Proteomics", "BRCA_proteomics_GO_BP_ORA")

## =============================================================================
## 6. PRAD PROTEOMICS — mass-spec (5,412 proteins, FULL) via ORA, custom background
##    NOTE: using NOMINAL p < 0.05 for significance (FDR < 0.05 leaves only 4 genes,
##    too few for ORA). Switch sig_col logic back to FDR if you'd rather keep that cut.
## =============================================================================
prad_dep <- read_delim("03_Results/Proteomics/prad_DEPs_all.tsv", delim = "\t", show_col_types = FALSE)

prad_dep <- prad_dep %>%
  filter(!is.na(GeneName)) %>%
  distinct(GeneName, .keep_all = TRUE)

prad_dep_map <- suppressWarnings(bitr(prad_dep$GeneName, "SYMBOL", "ENTREZID", org.Hs.eg.db))
prad_dep <- inner_join(prad_dep, prad_dep_map, by = c("GeneName" = "SYMBOL"))

prad_dep_universe <- unique(prad_dep$ENTREZID)
prad_dep_sig_fdr    <- unique(prad_dep$ENTREZID[prad_dep$FDR < 0.05])
prad_dep_sig_nominal <- unique(prad_dep$ENTREZID[prad_dep$p_value < 0.05])

cat("PRAD proteomics: ", length(prad_dep_sig_fdr), "sig at FDR<0.05 vs",
    length(prad_dep_sig_nominal), "at nominal p<0.05, of", length(prad_dep_universe), "background\n")

prad_dep_sig <- prad_dep_sig_nominal  # <-- change to prad_dep_sig_fdr if preferred

prad_dep_kegg <- enrichKEGG(gene = prad_dep_sig, universe = prad_dep_universe, organism = "hsa",
                             minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05)
save_enrich(prad_dep_kegg, "Proteomics", "PRAD_proteomics_KEGG_ORA")

prad_dep_go <- enrichGO(gene = prad_dep_sig, universe = prad_dep_universe, OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID", ont = "BP", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05)
save_enrich(prad_dep_go, "Proteomics", "PRAD_proteomics_GO_BP_ORA")

## =============================================================================
## Done — all outputs in 04_Pathway_Analysis/{Epigenomics,Transcriptomics,Proteomics}/
## =============================================================================
cat("\nStep 4.0 complete. Check console output above for BRCA RPPA mapping rate\n")
cat("and PRAD proteomics FDR-vs-nominal gene counts before interpreting results.\n")