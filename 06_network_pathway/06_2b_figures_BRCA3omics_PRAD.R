## =============================================================================
## STEP 6.2b — FIGURES for BRCA_3omics and PRAD_meth_rna (same treatment as
## BRCA_meth_rna, applied evenhandedly)
##
## IMPORTANT DIFFERENCE FROM BRCA_meth_rna: neither of these runs produced any
## FDR-significant (p.adjust <= 0.1) terms in ANY database, for either the
## FULL or STRICT gene set. So unlike BRCA_meth_rna's figure (which had a
## genuine two-theme, partly-significant story), these figures show the top
## NOMINAL terms only, explicitly labeled as such in the plot titles/legend —
## no color is used to imply significance that isn't there, since there is
## none to show. This is the same honest treatment, not a lesser one: we are
## not manufacturing a narrative these runs' statistics don't support.
##
## Requires: ggraph, tidygraph, patchwork (already installed for Step 6.2)
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
STEP6_DIR <- file.path(ROOT, "06_FinalResults", "Step6_Network")

RUNS <- c("BRCA_3omics", "PRAD_meth_rna")   # BRCA_meth_rna already done in Step 6.2


## ---- Helper: load a nominal top-N file, tag its source database ----
load_terms <- function(path, db_label) {
  if (!file.exists(path)) return(tibble())
  df <- read_csv(path, show_col_types = FALSE)
  if (nrow(df) == 0) return(tibble())
  df |>
    dplyr::transmute(term = Description, pvalue = pvalue, p.adjust = p.adjust,
                      Count = Count, db = db_label)
}

top_per_db <- function(df) {
  if (nrow(df) == 0) return(df)
  df |> dplyr::group_by(db) |> dplyr::slice_min(pvalue, n = 1) |> dplyr::ungroup()
}

make_panel <- function(data, panel_title) {
  if (nrow(data) == 0) {
    return(ggplot() + labs(title = panel_title, subtitle = "No enrichment terms available") + theme_void())
  }
  data <- data |> dplyr::mutate(term_short = stringr::str_trunc(term, 45),
                                  neg_log10_p = -log10(pvalue))
  ggplot(data, aes(x = neg_log10_p, y = reorder(term_short, neg_log10_p), size = Count)) +
    geom_point(color = "#7f7f7f") +
    labs(title = panel_title, x = expression(-log[10](p)), y = NULL, size = "Gene\ncount") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(size = 10, face = "bold"), legend.position = "right")
}


for (run_name in RUNS) {
  RUN_DIR <- file.path(STEP6_DIR, run_name)

  message(""); message("=== ", run_name, " ===")

  ## ---- Figure 1: FULL vs STRICT nominal top-term dotplot ----
  full_files <- list(
    Hallmark = file.path(RUN_DIR, sprintf("%s_FULL_Hallmark_top20_nominal.csv", run_name)),
    GO_BP    = file.path(RUN_DIR, sprintf("%s_FULL_GO_BP_top20_nominal.csv", run_name)),
    KEGG     = file.path(RUN_DIR, sprintf("%s_FULL_KEGG_top20_nominal.csv", run_name)),
    Reactome = file.path(RUN_DIR, sprintf("%s_FULL_Reactome_top20_nominal.csv", run_name))
  )
  full_terms <- bind_rows(lapply(names(full_files), function(db) load_terms(full_files[[db]], db)))
  full_top <- top_per_db(full_terms)

  strict_files <- list(
    Hallmark = file.path(RUN_DIR, sprintf("%s_STRICT_Hallmark_top20_nominal.csv", run_name)),
    GO_BP    = file.path(RUN_DIR, sprintf("%s_STRICT_GO_BP_top20_nominal.csv", run_name)),
    KEGG     = file.path(RUN_DIR, sprintf("%s_STRICT_KEGG_top20_nominal.csv", run_name)),
    Reactome = file.path(RUN_DIR, sprintf("%s_STRICT_Reactome_top20_nominal.csv", run_name))
  )
  strict_terms <- bind_rows(lapply(names(strict_files), function(db) load_terms(strict_files[[db]], db)))
  strict_top <- top_per_db(strict_terms)

  p1 <- make_panel(full_top,   sprintf("%s FULL set (nominal top terms, not FDR-significant)", run_name))
  p2 <- make_panel(strict_top, sprintf("%s STRICT set (nominal top terms, not FDR-significant)", run_name))

  fig1 <- p1 + p2 +
    patchwork::plot_annotation(
      title = sprintf("%s: exploratory nominal enrichment (no term reached p.adjust <= 0.1)", run_name),
      subtitle = "Top term per database (Hallmark, GO:BP, KEGG, Reactome) — reported for transparency, not as confirmed findings"
    )

  fig1_path <- file.path(RUN_DIR, sprintf("Figure1_%s_nominal_dotplot.png", run_name))
  ggsave(fig1_path, fig1, width = 12, height = 5, dpi = 300)
  message("Figure 1 saved: ", fig1_path)

  ## ---- Figure 2: STRING network for FULL set, colored by tier ----
  edges_path <- file.path(RUN_DIR, sprintf("%s_FULL_string_edges.csv", run_name))
  nodes_path <- file.path(RUN_DIR, sprintf("%s_FULL_node_attributes.csv", run_name))

  if (file.exists(edges_path) && file.exists(nodes_path)) {
    edges <- read_csv(edges_path, show_col_types = FALSE)
    nodes <- read_csv(nodes_path, show_col_types = FALSE)

    if (nrow(edges) == 0 || nrow(nodes) == 0) {
      message("!! ", run_name, ": empty network (0 nodes or edges) — skipping Figure 2")
    } else {
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
        labs(title = sprintf("%s — full consensus set STRING network", run_name),
             subtitle = "Node color = evidence tier | Node size = degree centrality",
             color = "Evidence tier", size = "Degree") +
        theme_void(base_size = 11) +
        theme(plot.title = element_text(face = "bold"), legend.position = "right")

      fig2_path <- file.path(RUN_DIR, sprintf("Figure2_%s_FULL_network.png", run_name))
      ggsave(fig2_path, fig2, width = 9, height = 7, dpi = 300)
      message("Figure 2 saved: ", fig2_path)
    }
  } else {
    message("!! ", run_name, ": missing edges/node attribute file — skipping Figure 2")
  }
}