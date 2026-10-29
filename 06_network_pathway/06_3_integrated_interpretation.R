## =============================================================================
## STEP 6.3 — INTEGRATED NETWORK-PATHWAY INTERPRETATION
##
## Cross-references the network analysis (hub genes) with the pathway
## enrichment results (which terms each gene belongs to), which the previous
## steps computed in parallel but never explicitly linked.
##
## Scope (proportionate to what each run's data supports):
##   - ALL THREE RUNS: hub-gene <-> pathway-term membership table (which
##     enriched/nominal terms each network hub gene belongs to)
##   - BRCA_meth_rna ONLY (the one run with real FDR-significant enrichment):
##     a full gene-concept network (cnetplot) linking genes to their
##     significant terms, with hub genes indicated
##
## Namespacing: dplyr/igraph calls explicitly prefixed throughout, per the
## collisions hit earlier in this pipeline (AnnotationDbi/S4Vectors shadow
## dplyr::select/filter and igraph::simplify once loaded).
## =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(tibble)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(msigdbr)
  library(enrichplot)
  library(ggplot2)
})
library(dplyr)  # loaded but always called with dplyr:: prefix below

ROOT      <- "C:/Users/USER/Documents/MyThesis_MultiOmics"
STEP6_DIR <- file.path(ROOT, "06_FinalResults", "Step6_Network")
RUNS      <- c("BRCA_3omics", "BRCA_meth_rna", "PRAD_meth_rna")

report <- c()
add <- function(...) report <<- c(report, sprintf(...))


## =============================================================================
## PART A — Hub-gene <-> pathway membership table (all three runs)
## =============================================================================

## Split a "/"-joined geneID string into a vector, mapping Entrez -> symbol
## when needed (KEGG results are Entrez-coded; GO/Reactome/Hallmark are
## already gene symbols in our pipeline since readable=TRUE / msigdbr symbols
## were used).
genes_from_geneID <- function(geneID_str, is_entrez = FALSE) {
  ids <- stringr::str_split(geneID_str, "/")[[1]]
  if (!is_entrez) return(ids)
  suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, keys = ids, column = "SYMBOL",
                                          keytype = "ENTREZID", multiVals = "first")) |>
    unname() |> stats::na.omit() |> as.character()
}

## Build a long (term, gene, db, pvalue) table from a top-N nominal CSV
term_gene_table <- function(path, db_label, is_entrez = FALSE) {
  if (!file.exists(path)) return(tibble())
  df <- read_csv(path, show_col_types = FALSE)
  if (nrow(df) == 0 || !"geneID" %in% names(df)) return(tibble())
  purrr_rows <- lapply(seq_len(nrow(df)), function(i) {
    genes <- genes_from_geneID(df$geneID[i], is_entrez = is_entrez)
    if (length(genes) == 0) return(tibble())
    tibble(term = df$Description[i], gene = genes, db = db_label,
           pvalue = df$pvalue[i], p.adjust = df$p.adjust[i])
  })
  dplyr::bind_rows(purrr_rows)
}

for (run_name in RUNS) {
  outdir <- file.path(STEP6_DIR, run_name)
  add(""); add("=====================================================================")
  add("RUN: %s", run_name)
  add("=====================================================================")

  centrality_path <- file.path(outdir, sprintf("%s_FULL_network_centrality.csv", run_name))
  if (!file.exists(centrality_path)) { add("  !! no network centrality file — skipping"); next }
  centrality <- read_csv(centrality_path, show_col_types = FALSE)

  ## Gather term-gene tables from all four databases (FULL set)
  tg <- dplyr::bind_rows(
    term_gene_table(file.path(outdir, sprintf("%s_FULL_Hallmark_top20_nominal.csv", run_name)), "Hallmark"),
    term_gene_table(file.path(outdir, sprintf("%s_FULL_GO_BP_top20_nominal.csv", run_name)), "GO_BP"),
    term_gene_table(file.path(outdir, sprintf("%s_FULL_KEGG_top20_nominal.csv", run_name)), "KEGG", is_entrez = TRUE),
    term_gene_table(file.path(outdir, sprintf("%s_FULL_Reactome_top20_nominal.csv", run_name)), "Reactome")
  )

  if (nrow(tg) == 0) {
    add("  !! no gene-term membership data available (no enrichment terms found)")
    next
  }

  ## For each network gene (not just top hubs — small networks, so show all),
  ## list which terms it belongs to
  gene_term_summary <- tg |>
    dplyr::group_by(gene) |>
    dplyr::summarise(
      n_terms = dplyr::n(),
      terms = paste(sprintf("%s: %s (p=%.3g)", db, term, pvalue), collapse = " | "),
      any_fdr_sig = any(p.adjust <= 0.1),
      .groups = "drop"
    )

  hub_pathway <- centrality |>
    dplyr::left_join(gene_term_summary, by = "gene") |>
    dplyr::mutate(n_terms = tidyr::replace_na(n_terms, 0),
                   any_fdr_sig = tidyr::replace_na(any_fdr_sig, FALSE)) |>
    dplyr::arrange(dplyr::desc(degree))

  write_csv(hub_pathway, file.path(outdir, sprintf("%s_hub_pathway_overlap.csv", run_name)))

  top_hub <- hub_pathway$gene[1]
  hub_in_pathway <- hub_pathway$n_terms[1] > 0
  add("  Top hub gene: %s (degree %d) — %s", top_hub, hub_pathway$degree[1],
      ifelse(hub_in_pathway, "IS a member of an enriched/nominal term", "is NOT in any captured enriched term"))
  n_hubs_with_terms <- sum(hub_pathway$n_terms > 0)
  add("  %d / %d network genes appear in at least one captured pathway term",
      n_hubs_with_terms, nrow(hub_pathway))
  add("  Table written -> %s", file.path(outdir, sprintf("%s_hub_pathway_overlap.csv", run_name)))
}


## =============================================================================
## PART B — Full gene-concept network (cnetplot) for BRCA_meth_rna only
## =============================================================================

add(""); add("=====================================================================")
add("PART B: BRCA_meth_rna gene-concept network (significant terms only)")
add("=====================================================================")

run_name <- "BRCA_meth_rna"
outdir <- file.path(STEP6_DIR, run_name)
cg_path <- file.path(outdir, "consensus_genes.csv")

if (file.exists(cg_path)) {
  cg <- read_csv(cg_path, show_col_types = FALSE)
  all_genes <- unique(cg$gene)
  entrez <- suppressMessages(AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = all_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"
  )) |> stats::na.omit()

  ## Re-run the two significant enrichments as real enrichResult objects
  ## (needed for cnetplot; only the data frames were saved previously)
  go_bp <- tryCatch(clusterProfiler::enrichGO(entrez, org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                                               pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE),
                     error = function(e) NULL)
  hallmark_sets <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", category = "H"), error = function(e) NULL)

  if (!is.null(go_bp)) {
    go_df <- as.data.frame(go_bp)
    sig_go <- go_df |> dplyr::filter(p.adjust <= 0.1)
    add("  GO:BP significant terms for cnetplot: %d", nrow(sig_go))

    if (nrow(sig_go) > 0) {
      ## Restrict the enrichResult to just the significant term(s) for a focused plot
      go_bp_sig <- go_bp
      go_bp_sig@result <- sig_go

      fig <- enrichplot::cnetplot(go_bp_sig, showCategory = nrow(sig_go), node_label = "all") +
        ggtitle("BRCA methylation+RNA — gene-concept network (GO:BP, FDR-significant terms)")

      fig_path <- file.path(outdir, "Figure3_BRCA_methrna_cnetplot_GO.png")
      ggsave(fig_path, fig, width = 9, height = 7, dpi = 300)
      add("  cnetplot saved -> %s", fig_path)
    }
  } else {
    add("  !! GO:BP re-run failed — skipping cnetplot")
  }
} else {
  add("  !! consensus_genes.csv not found for BRCA_meth_rna — skipping Part B")
}

REPORT_OUT <- file.path(STEP6_DIR, "Step6_3_integrated_interpretation_report.txt")
writeLines(report, REPORT_OUT)
cat(paste(report, collapse = "\n"))
cat("\n\nFull report saved to:\n", REPORT_OUT, "\n")