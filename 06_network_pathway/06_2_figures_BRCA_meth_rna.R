## =============================================================================
## STEP 6.2 — FIGURES: BRCA_meth_rna immune vs ECM signature comparison,
##            and STRING network colored by consensus tier
##
## Requires (one-time install if not already present):
##   install.packages(c("ggraph", "tidygraph", "patchwork"))
##
## Produces two figures in 06_FinalResults/Step6_Network/BRCA_meth_rna/:
##   Figure1_Immune_vs_ECM_dotplot.png
##       Left panel  = FULL consensus set (immune/T-cell theme) top terms
##       Right panel = STRICT (core+moderate) set (ECM theme) top terms
##       Bars/dots colored by whether p.adjust <= 0.1 (FDR-significant) or not
##   Figure2_BRCA_methrna_FULL_network.png
##       STRING network for the FULL consensus set, nodes colored by tier,
##       sized by degree, labeled with gene symbol
## =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggraph)
  library(tidygraph)
  library(patchwork)
})

ROOT      <- "C:/Users/USER/Documents/MyThesis_MultiOmics"
RUN_DIR   <- file.path(ROOT, "06_FinalResults", "Step6_Network", "BRCA_meth_rna")

## =============================================================================
## FIGURE 1 — Immune (FULL) vs ECM (STRICT) top-term dotplot
## =============================================================================

## ---- Helper: load a nominal top-N file, tag its source database ----
load_terms <- function(path, db_label) {
  if (!file.exists(path)) return(tibble())
  df <- read_csv(path, show_col_types = FALSE)
  if (nrow(df) == 0) return(tibble())
  df |>
    dplyr::transmute(
      term = Description,
      pvalue = pvalue,
      p.adjust = p.adjust,
      Count = Count,
      db = db_label
    )
}

## Curate the specific terms that told the coherent story in the console
## output (rather than dumping all 20 nominal terms per file, which would
## bury the signal) — edit the `keep_terms` vectors below if you'd rather
## show a different subset once you've reviewed the full nominal CSVs.

full_files <- list(
  Hallmark = file.path(RUN_DIR, "BRCA_meth_rna_FULL_Hallmark_top20_nominal.csv"),
  GO_BP    = file.path(RUN_DIR, "BRCA_meth_rna_FULL_GO_BP_top20_nominal.csv"),
  KEGG     = file.path(RUN_DIR, "BRCA_meth_rna_FULL_KEGG_top20_nominal.csv"),
  Reactome = file.path(RUN_DIR, "BRCA_meth_rna_FULL_Reactome_top20_nominal.csv")
)
full_terms <- bind_rows(lapply(names(full_files), function(db) load_terms(full_files[[db]], db)))

strict_files <- list(
  Hallmark = file.path(RUN_DIR, "BRCA_meth_rna_STRICT_Hallmark_top20_nominal.csv"),
  GO_BP    = file.path(RUN_DIR, "BRCA_meth_rna_STRICT_GO_BP_top20_nominal.csv"),
  KEGG     = file.path(RUN_DIR, "BRCA_meth_rna_STRICT_KEGG_top20_nominal.csv"),
  Reactome = file.path(RUN_DIR, "BRCA_meth_rna_STRICT_Reactome_top20_nominal.csv")
)
strict_terms <- bind_rows(lapply(names(strict_files), function(db) load_terms(strict_files[[db]], db)))

## Take the single top term per database for each panel (this reproduces
## exactly the 4 terms highlighted in the console summary for each set).
top_per_db <- function(df) {
  if (nrow(df) == 0) return(df)
  df |> dplyr::group_by(db) |> dplyr::slice_min(pvalue, n = 1) |> dplyr::ungroup()
}

full_top   <- top_per_db(full_terms)   |> dplyr::mutate(panel = "FULL consensus set (immune/T-cell theme)")
strict_top <- top_per_db(strict_terms) |> dplyr::mutate(panel = "STRICT core+moderate set (ECM theme)")

plot_df <- dplyr::bind_rows(full_top, strict_top) |>
  dplyr::mutate(
    term_short = stringr::str_trunc(term, 45),
    sig = ifelse(p.adjust <= 0.1, "FDR p.adjust <= 0.1", "Nominal only"),
    neg_log10_p = -log10(pvalue)
  )

make_panel <- function(data, panel_title) {
  ggplot(data, aes(x = neg_log10_p, y = reorder(term_short, neg_log10_p), size = Count, color = sig)) +
    geom_point() +
    scale_color_manual(values = c("FDR p.adjust <= 0.1" = "#d62728", "Nominal only" = "#7f7f7f")) +
    labs(title = panel_title, x = expression(-log[10](p)), y = NULL, size = "Gene\ncount", color = "Significance") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(size = 10, face = "bold"), legend.position = "right")
}

p1 <- make_panel(dplyr::filter(plot_df, panel == "FULL consensus set (immune/T-cell theme)"),
                  "FULL set: Immune/T-cell signature")
p2 <- make_panel(dplyr::filter(plot_df, panel == "STRICT core+moderate set (ECM theme)"),
                  "STRICT set: ECM signature")

fig1 <- p1 + p2 + patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title = "BRCA methylation+RNA integration: two distinct signatures by evidence tier",
    subtitle = "Top term per database (Hallmark, GO:BP, KEGG, Reactome); color = FDR-significance"
  )

ggsave(file.path(RUN_DIR, "Figure1_Immune_vs_ECM_dotplot.png"), fig1, width = 12, height = 5, dpi = 300)
message("Figure 1 saved: ", file.path(RUN_DIR, "Figure1_Immune_vs_ECM_dotplot.png"))


## =============================================================================
## FIGURE 2 — STRING network for BRCA_meth_rna FULL, colored by tier
## =============================================================================

edges_path <- file.path(RUN_DIR, "BRCA_meth_rna_FULL_string_edges.csv")
nodes_path <- file.path(RUN_DIR, "BRCA_meth_rna_FULL_node_attributes.csv")

if (file.exists(edges_path) && file.exists(nodes_path)) {
  edges <- read_csv(edges_path, show_col_types = FALSE)
  nodes <- read_csv(nodes_path, show_col_types = FALSE)

  ## Simplify tier labels for legend readability
  nodes <- nodes |>
    dplyr::mutate(tier_short = dplyr::case_when(
      stringr::str_detect(tier, "core") ~ "Core (DIABLO+DE+MOFA)",
      stringr::str_detect(tier, "moderate") ~ "Moderate (DIABLO+1)",
      TRUE ~ "Extended (DIABLO only)"
    ))

  graph_tbl <- tidygraph::tbl_graph(
    nodes = nodes |> dplyr::rename(name = gene),
    edges = edges |> dplyr::select(gene_from, gene_to) |> dplyr::rename(from = gene_from, to = gene_to),
    node_key = "name", directed = FALSE
  )

  fig2 <- ggraph::ggraph(graph_tbl, layout = "fr") +
    ggraph::geom_edge_link(alpha = 0.3, color = "grey60") +
    ggraph::geom_node_point(aes(size = degree, color = tier_short)) +
    ggraph::geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    scale_color_manual(values = c(
      "Core (DIABLO+DE+MOFA)" = "#d62728",
      "Moderate (DIABLO+1)"   = "#ff7f0e",
      "Extended (DIABLO only)" = "#1f77b4"
    )) +
    labs(title = "BRCA methylation+RNA — full consensus set STRING network",
         subtitle = "Node color = evidence tier | Node size = degree centrality",
         color = "Evidence tier", size = "Degree") +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")

  ggsave(file.path(RUN_DIR, "Figure2_BRCA_methrna_FULL_network.png"), fig2, width = 9, height = 7, dpi = 300)
  message("Figure 2 saved: ", file.path(RUN_DIR, "Figure2_BRCA_methrna_FULL_network.png"))
} else {
  message("!! Missing edges or node attribute file — check paths:\n  ", edges_path, "\n  ", nodes_path)
}