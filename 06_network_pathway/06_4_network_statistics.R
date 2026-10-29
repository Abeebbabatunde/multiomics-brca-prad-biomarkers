## =============================================================================
## STEP 6.4 — FULL NETWORK STATISTICS + CONSENSUS HUB RANKING
##
## Adds two things missing from the per-gene centrality tables already
## produced (Step 6.1):
##   1. Whole-network descriptive statistics: density, average degree,
##      clustering coefficient (transitivity), number of connected
##      components, diameter, average path length — computed on the FULL
##      consensus-set network for each run.
##   2. Multi-measure consensus hub ranking: rather than trusting degree
##      alone, genes are ranked separately by degree, betweenness, closeness,
##      and eigenvector centrality; a consensus rank (mean of the four ranks)
##      identifies genes that are robustly central across measures, not just
##      by one metric.
##
## Namespacing: igraph:: prefixed throughout (clusterProfiler's S4 'simplify'
## generic shadows igraph::simplify once both are loaded — a collision hit
## earlier in this pipeline).
## =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(tibble)
})
library(dplyr)   # always called with dplyr:: prefix below
library(igraph)  # always called with igraph:: prefix below

ROOT      <- "C:/Users/USER/Documents/MyThesis_MultiOmics"
STEP6_DIR <- file.path(ROOT, "06_FinalResults", "Step6_Network")
RUNS      <- c("BRCA_3omics", "BRCA_meth_rna", "PRAD_meth_rna")

report <- c()
add <- function(...) report <<- c(report, sprintf(...))

network_summary_rows <- list()

for (run_name in RUNS) {
  outdir <- file.path(STEP6_DIR, run_name)
  edges_path <- file.path(outdir, sprintf("%s_FULL_string_edges.csv", run_name))
  centrality_path <- file.path(outdir, sprintf("%s_FULL_network_centrality.csv", run_name))

  add(""); add("=====================================================================")
  add("RUN: %s", run_name)
  add("=====================================================================")

  if (!file.exists(edges_path)) { add("  !! edges file not found — skipping"); next }
  edges <- read_csv(edges_path, show_col_types = FALSE)

  if (nrow(edges) == 0) { add("  !! network has 0 edges — no statistics to compute"); next }

  g <- igraph::simplify(igraph::graph_from_data_frame(
    edges |> dplyr::select(gene_from, gene_to), directed = FALSE
  ))

  n_nodes <- igraph::vcount(g)
  n_edges <- igraph::ecount(g)
  density <- igraph::edge_density(g)
  avg_degree <- mean(igraph::degree(g))
  clustering_coef <- igraph::transitivity(g, type = "average", isolates = "zero")
  comp <- igraph::components(g)
  n_components <- comp$no
  largest_comp_size <- max(comp$csize)

  ## Diameter / average path length computed over the largest connected
  ## component only (unconnected=TRUE ignores unreachable pairs, which
  ## otherwise return Inf for a disconnected graph)
  diam <- tryCatch(igraph::diameter(g, unconnected = TRUE), error = function(e) NA_real_)
  avg_path_len <- tryCatch(igraph::mean_distance(g, unconnected = TRUE), error = function(e) NA_real_)

  add("  Nodes: %d | Edges: %d | Density: %.4f", n_nodes, n_edges, density)
  add("  Average degree: %.2f | Clustering coefficient (avg. transitivity): %.4f", avg_degree, clustering_coef)
  add("  Connected components: %d | Largest component size: %d nodes", n_components, largest_comp_size)
  add("  Diameter (largest component): %s | Average path length: %s",
      ifelse(is.na(diam), "NA", diam), ifelse(is.na(avg_path_len), "NA", round(avg_path_len, 3)))

  network_summary_rows[[run_name]] <- tibble(
    run = run_name, n_nodes = n_nodes, n_edges = n_edges, density = density,
    avg_degree = avg_degree, clustering_coefficient = clustering_coef,
    n_components = n_components, largest_component_size = largest_comp_size,
    diameter = diam, avg_path_length = avg_path_len
  )

  ## ---- Consensus hub ranking ----
  if (!file.exists(centrality_path)) { add("  !! centrality file not found — skipping hub ranking"); next }
  centrality <- read_csv(centrality_path, show_col_types = FALSE)

  if (nrow(centrality) < 2) {
    add("  !! too few nodes (%d) for a meaningful consensus ranking", nrow(centrality))
    next
  }

  ranked <- centrality |>
    dplyr::mutate(
      rank_degree      = rank(-degree, ties.method = "average"),
      rank_betweenness  = rank(-betweenness, ties.method = "average"),
      rank_closeness    = rank(-closeness, ties.method = "average"),
      rank_eigenvector  = rank(-eigenvector, ties.method = "average"),
      consensus_rank    = (rank_degree + rank_betweenness + rank_closeness + rank_eigenvector) / 4,
      top3_count = (rank_degree <= 3) + (rank_betweenness <= 3) + (rank_closeness <= 3) + (rank_eigenvector <= 3)
    ) |>
    dplyr::arrange(consensus_rank)

  write_csv(ranked, file.path(outdir, sprintf("%s_FULL_consensus_hub_ranking.csv", run_name)))

  top_consensus <- ranked$gene[1]
  top_degree_only <- centrality$gene[which.max(centrality$degree)]
  agree <- top_consensus == top_degree_only

  add("  Top hub by DEGREE alone: %s", top_degree_only)
  add("  Top hub by CONSENSUS (avg rank across 4 measures): %s%s",
      top_consensus, ifelse(agree, " (agrees with degree-only pick)", " (DIFFERS from degree-only pick)"))

  robust_hubs <- ranked |> dplyr::filter(top3_count >= 3) |> dplyr::pull(gene)
  if (length(robust_hubs) > 0) {
    add("  Genes ranking in the top 3 on >=3 of 4 centrality measures (most defensible hubs): %s",
        paste(robust_hubs, collapse = ", "))
  } else {
    add("  No gene ranks in the top 3 on 3+ of the 4 measures simultaneously — no strongly consensus hub")
  }
}

## ---- Combined cross-run network summary table ----
if (length(network_summary_rows) > 0) {
  combined <- dplyr::bind_rows(network_summary_rows)
  write_csv(combined, file.path(STEP6_DIR, "cross_run_network_statistics.csv"))
  add(""); add("Combined cross-run network statistics table -> %s",
      file.path(STEP6_DIR, "cross_run_network_statistics.csv"))
}

REPORT_OUT <- file.path(STEP6_DIR, "Step6_4_network_statistics_report.txt")
writeLines(report, REPORT_OUT)
cat(paste(report, collapse = "\n"))
cat("\n\nFull report saved to:\n", REPORT_OUT, "\n")