## =============================================================================
## STEP 6.1 — INTEGRATED NETWORK + PATHWAY ANALYSIS (real data) — CORRECTED
##
## Fixes vs previous version:
##   - options(timeout = 600) to avoid STRINGdb download timeouts on large
##     alias/link files
##   - All igraph calls explicitly namespaced (igraph::simplify, igraph::degree,
##     igraph::betweenness, igraph::closeness, igraph::eigen_centrality,
##     igraph::vcount, igraph::ecount, igraph::V, igraph::graph_from_data_frame)
##     because clusterProfiler/BiocGenerics defines its own S4 'simplify'
##     generic that silently shadows igraph::simplify once both packages are
##     loaded, causing "unable to find an inherited method for function
##     'simplify' for signature x = 'igraph'".
##   - dplyr verbs explicitly namespaced (dplyr::rename/filter/select/...)
##     because AnnotationDbi::select masks dplyr::select.
##   - get_interactions() already returns from_name / to_name gene symbols, so
##     they are renamed directly (no id2gene lookup; the old id2gene[from]
##     failed because 'from' resolved to igraph::from).
##
## Input: the three consensus_genes.csv files built in Step 6.0a/6.0c.
## =============================================================================

options(timeout = 600)

library(dplyr)
library(readr)
library(tibble)
library(ggplot2)
library(STRINGdb)
library(igraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

ROOT      <- "C:/Users/USER/Documents/MyThesis_MultiOmics"
STEP6_DIR <- file.path(ROOT, "06_FinalResults", "Step6_Network")

RUNS <- c("BRCA_3omics", "BRCA_meth_rna", "PRAD_meth_rna")

STRING_SCORE_THRESHOLD <- 400   # medium confidence; use 700 for high-confidence/publication

report <- c()
add <- function(...) report <<- c(report, sprintf(...))


## ---- Core function: STRING network + hubs -------------------------------------

build_network_and_hubs <- function(genes, outdir, label) {
  string_db <- STRINGdb$new(version = "12.0", species = 9606,
                             score_threshold = STRING_SCORE_THRESHOLD,
                             input_directory = outdir)
  gene_df <- data.frame(gene = genes, stringsAsFactors = FALSE)
  mapped  <- string_db$map(gene_df, "gene", removeUnmappedRows = TRUE)
  add("  [%s] mapped %d / %d genes to STRING", label, nrow(mapped), length(genes))

  if (nrow(mapped) < 3) {
    add("  !! too few mapped genes for a meaningful network — skipping network step")
    return(NULL)
  }

  edges <- string_db$get_interactions(mapped$STRING_id)

  if (is.null(edges) || nrow(edges) == 0) {
    add("  !! STRING returned 0 interactions at threshold %d — network is empty for this set",
        STRING_SCORE_THRESHOLD)
    return(list(mapped = mapped, edges = edges, graph = NULL, centrality = NULL))
  }

  ## get_interactions() already provides readable gene symbols directly —
  ## no manual STRING_id -> gene mapping needed (and safer: avoids a silent
  ## fallback to an unrelated same-named function if a column were missing).
  edges <- edges |>
    dplyr::rename(gene_from = from_name, gene_to = to_name) |>
    dplyr::filter(!is.na(gene_from), !is.na(gene_to))

  write_csv(edges, file.path(outdir, sprintf("%s_string_edges.csv", label)))

  g <- igraph::simplify(igraph::graph_from_data_frame(edges[, c("gene_from", "gene_to")], directed = FALSE))

  centrality <- tibble(
    gene = igraph::V(g)$name,
    degree = igraph::degree(g),
    betweenness = igraph::betweenness(g, normalized = TRUE),
    closeness = igraph::closeness(g, normalized = TRUE),
    eigenvector = igraph::eigen_centrality(g)$vector
  ) |> dplyr::arrange(dplyr::desc(degree))

  write_csv(centrality, file.path(outdir, sprintf("%s_network_centrality.csv", label)))
  add("  [%s] network: %d nodes, %d edges. Top hub: %s (degree %d)",
      label, igraph::vcount(g), igraph::ecount(g),
      centrality$gene[1], centrality$degree[1])

  list(mapped = mapped, edges = edges, graph = g, centrality = centrality)
}


## ---- Core function: enrichment suite ------------------------------------------

run_enrichment <- function(genes, outdir, label) {
  entrez <- suppressMessages(AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"
  )) |> na.omit()

  if (length(entrez) < 3) {
    add("  !! too few mapped Entrez IDs (%d) for %s — skipping enrichment", length(entrez), label)
    return(NULL)
  }

  go_bp    <- tryCatch(enrichGO(entrez, org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                                 pvalueCutoff = 0.1, qvalueCutoff = 0.2, readable = TRUE),
                        error = function(e) NULL)
  kegg     <- tryCatch(enrichKEGG(entrez, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.1),
                        error = function(e) NULL)
  reactome <- tryCatch(enrichPathway(entrez, organism = "human", pAdjustMethod = "BH",
                                      pvalueCutoff = 0.1, readable = TRUE),
                        error = function(e) NULL)

  n_go  <- if (!is.null(go_bp))    nrow(as.data.frame(go_bp))    else 0
  n_kg  <- if (!is.null(kegg))     nrow(as.data.frame(kegg))     else 0
  n_rc  <- if (!is.null(reactome)) nrow(as.data.frame(reactome)) else 0
  add("  [%s] enrichment hits (p.adj<=0.1): GO:BP=%d, KEGG=%d, Reactome=%d", label, n_go, n_kg, n_rc)

  if (n_go > 0) write_csv(as.data.frame(go_bp),    file.path(outdir, sprintf("%s_GO_BP.csv", label)))
  if (n_kg > 0) write_csv(as.data.frame(kegg),     file.path(outdir, sprintf("%s_KEGG.csv", label)))
  if (n_rc > 0) write_csv(as.data.frame(reactome), file.path(outdir, sprintf("%s_Reactome.csv", label)))

  if (n_go > 0) {
    ggsave(file.path(outdir, sprintf("%s_GO_BP_dotplot.png", label)),
           dotplot(go_bp, showCategory = 15) + ggtitle(label), width = 8, height = 7, dpi = 300)
  }

  list(go_bp = go_bp, kegg = kegg, reactome = reactome)
}


## ---- Cytoscape export ----------------------------------------------------------

export_cytoscape <- function(edges, centrality, tier_map, outdir, label) {
  if (is.null(edges) || nrow(edges) == 0) return(invisible(NULL))
  sif <- edges |> dplyr::mutate(interaction = "pp") |> dplyr::select(gene_from, interaction, gene_to)
  write_delim(sif, file.path(outdir, sprintf("%s_network.sif", label)), delim = "\t", col_names = FALSE)

  node_attrs <- centrality |> dplyr::left_join(tier_map, by = "gene")
  write_csv(node_attrs, file.path(outdir, sprintf("%s_node_attributes.csv", label)))
}


## ---- Main loop per run -----------------------------------------------------------

enrichment_by_run <- list()

for (run_name in RUNS) {
  outdir <- file.path(STEP6_DIR, run_name)
  cg_path <- file.path(outdir, "consensus_genes.csv")

  add(""); add("=====================================================================")
  add("RUN: %s", run_name)
  add("=====================================================================")

  if (!file.exists(cg_path)) { add("  !! consensus_genes.csv not found — skipping"); next }
  cg <- read_csv(cg_path, show_col_types = FALSE)

  all_genes <- unique(cg$gene)
  strict_genes <- cg |> dplyr::filter(tier != "extended (DIABLO only)") |> dplyr::pull(gene) |> unique()

  tier_map <- cg |> dplyr::group_by(gene) |>
    dplyr::summarise(tier = paste(unique(tier), collapse = "; "), n_support = max(n_support), .groups = "drop")

  add("Full consensus set: %d genes | Core+moderate (strict) subset: %d genes",
      length(all_genes), length(strict_genes))

  ## --- Primary: full consensus set ---
  add(""); add("-- PRIMARY (full consensus set) --")
  net_full <- build_network_and_hubs(all_genes, outdir, paste0(run_name, "_FULL"))
  enr_full <- run_enrichment(all_genes, outdir, paste0(run_name, "_FULL"))
  if (!is.null(net_full) && !is.null(net_full$centrality)) {
    export_cytoscape(net_full$edges, net_full$centrality, tier_map, outdir, paste0(run_name, "_FULL"))
  }

  ## --- Sensitivity: core + moderate only ---
  add(""); add("-- SENSITIVITY (core + moderate only, DIABLO-only genes excluded) --")
  if (length(strict_genes) >= 3) {
    net_strict <- build_network_and_hubs(strict_genes, outdir, paste0(run_name, "_STRICT"))
    enr_strict <- run_enrichment(strict_genes, outdir, paste0(run_name, "_STRICT"))
  } else {
    add("  !! too few strict-tier genes (%d) for a separate sensitivity network/enrichment", length(strict_genes))
    enr_strict <- NULL
  }

  enrichment_by_run[[run_name]] <- list(full = enr_full, strict = enr_strict)
}


## ---- Cross-run GO:BP comparison (full consensus set) ------------------------------

add(""); add("=====================================================================")
add("CROSS-RUN GO:BP COMPARISON (full consensus sets)")
add("=====================================================================")

go_terms <- lapply(enrichment_by_run, function(x) {
  if (!is.null(x$full) && !is.null(x$full$go_bp)) as.data.frame(x$full$go_bp)$Description else character()
})

if (length(go_terms) >= 2 && all(sapply(go_terms, length) > 0)) {
  all_terms <- unique(unlist(go_terms))
  comparison <- tibble(term = all_terms)
  for (rn in names(go_terms)) comparison[[rn]] <- comparison$term %in% go_terms[[rn]]
  write_csv(comparison, file.path(STEP6_DIR, "cross_run_GO_BP_comparison.csv"))
  shared_all <- comparison |> dplyr::filter(dplyr::if_all(dplyr::all_of(names(go_terms)), ~ .x)) |> nrow()
  add("Terms shared across ALL runs with enrichment: %d", shared_all)
  add("Full comparison table -> %s", file.path(STEP6_DIR, "cross_run_GO_BP_comparison.csv"))
} else {
  add("  Not enough runs with GO:BP hits to compare (check enrichment hit counts above).")
}

REPORT_OUT <- file.path(STEP6_DIR, "Step6_1_network_pathway_report.txt")
writeLines(report, REPORT_OUT)
cat(paste(report, collapse = "\n"))
cat("\n\nFull report saved to:\n", REPORT_OUT, "\n")
