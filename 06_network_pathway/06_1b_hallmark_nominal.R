## =============================================================================
## STEP 6.1b — HALLMARK ENRICHMENT + NOMINAL TOP TERMS
##
## Adds two things the small consensus panels (28-48 genes) need:
##   1. MSigDB Hallmark enrichment (50 curated cancer-relevant gene sets) via
##      msigdbr + clusterProfiler::enricher() — far less multiple-testing
##      burden than full GO (~12,000+ terms), better suited to small panels.
##      Uses gene SYMBOLS directly (msigdbr provides them), avoiding the
##      Entrez-mapping loss seen in the GO/KEGG/Reactome step.
##   2. Nominal (uncorrected) top-N terms for GO:BP/KEGG/Reactome, exported
##      regardless of BH-adjusted significance, clearly labeled as exploratory
##      — standard practice for small-panel work, as long as it's disclosed
##      rather than presented as FDR-significant.
##
## Namespacing note (lessons from this session's crashes): dplyr verbs
## (select/mutate/filter/arrange/desc/rename/pull/group_by/summarise/
## left_join/bind_rows) and igraph functions are called with explicit
## dplyr:: / igraph:: prefixes throughout, because AnnotationDbi/S4Vectors/
## clusterProfiler define same-named S4 generics (select, simplify, from, to)
## that silently shadow the intended function once loaded — this caused three
## separate crashes earlier in this pipeline. Do not remove the prefixes.
## =============================================================================

options(timeout = 600)

suppressPackageStartupMessages({
  library(readr)
  library(tibble)
  library(msigdbr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
})
# dplyr loaded but NEVER called unqualified below — see note above
library(dplyr)

ROOT      <- "C:/Users/USER/Documents/MyThesis_MultiOmics"
STEP6_DIR <- file.path(ROOT, "06_FinalResults", "Step6_Network")
RUNS      <- c("BRCA_3omics", "BRCA_meth_rna", "PRAD_meth_rna")
N_NOMINAL <- 20   # top-N nominal terms to export regardless of adjusted significance

report <- c()
add <- function(...) report <<- c(report, sprintf(...))


## ---- Load Hallmark gene sets once (shared across all runs) -------------------

hallmark_sets <- tryCatch({
  msigdbr::msigdbr(species = "Homo sapiens", category = "H")
}, error = function(e) {
  add("!! msigdbr failed to load Hallmark sets: %s", conditionMessage(e))
  NULL
})

if (!is.null(hallmark_sets)) {
  # clusterProfiler::enricher() TERM2GENE expects two columns: term, gene
  gene_col <- if ("gene_symbol" %in% names(hallmark_sets)) "gene_symbol" else "human_gene_symbol"
  hallmark_t2g <- hallmark_sets |>
    dplyr::select(gs_name, dplyr::all_of(gene_col)) |>
    dplyr::rename(term = gs_name, gene = dplyr::all_of(gene_col)) |>
    dplyr::distinct()
  add("Hallmark gene sets loaded: %d gene sets, %d gene-term pairs",
      length(unique(hallmark_t2g$term)), nrow(hallmark_t2g))
}


## ---- Hallmark enrichment (gene symbols, no Entrez mapping needed) -------------

run_hallmark <- function(genes, outdir, label) {
  if (is.null(hallmark_sets)) { add("  !! Hallmark sets unavailable — skipping %s", label); return(NULL) }

  res <- tryCatch(
    clusterProfiler::enricher(
      gene = genes, TERM2GENE = hallmark_t2g,
      pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1  # capture everything; filter after
    ),
    error = function(e) { add("  !! enricher() failed for %s: %s", label, conditionMessage(e)); NULL }
  )

  if (is.null(res)) return(NULL)
  res_df <- as.data.frame(res) |> dplyr::arrange(pvalue)

  if (nrow(res_df) == 0) { add("  [%s] Hallmark: 0 gene sets returned any overlap", label); return(NULL) }

  n_sig <- sum(res_df$p.adjust <= 0.1, na.rm = TRUE)
  add("  [%s] Hallmark: %d gene sets tested, %d with p.adjust<=0.1, top nominal: %s (p=%.4g)",
      label, nrow(res_df), n_sig, res_df$Description[1], res_df$pvalue[1])

  write_csv(res_df, file.path(outdir, sprintf("%s_Hallmark_full.csv", label)))
  write_csv(utils::head(res_df, N_NOMINAL), file.path(outdir, sprintf("%s_Hallmark_top%d_nominal.csv", label, N_NOMINAL)))
  res_df
}


## ---- GO/KEGG/Reactome re-run with nominal top-N export -------------------------

run_enrichment_with_nominal <- function(genes, outdir, label) {
  entrez <- suppressMessages(AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"
  ))
  entrez <- entrez[!is.na(entrez)]

  if (length(entrez) < 3) {
    add("  !! too few mapped Entrez IDs (%d) for %s — skipping", length(entrez), label)
    return(invisible(NULL))
  }

  run_one <- function(fun, name, ...) {
    res <- tryCatch(fun(entrez, ..., pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1),
                     error = function(e) { add("  !! %s failed for %s: %s", name, label, conditionMessage(e)); NULL })
    if (is.null(res)) return(invisible(NULL))
    df <- as.data.frame(res)
    if (nrow(df) == 0) { add("  [%s] %s: 0 terms returned", label, name); return(invisible(NULL)) }
    df <- df |> dplyr::arrange(pvalue)
    n_sig <- sum(df$p.adjust <= 0.1, na.rm = TRUE)
    add("  [%s] %s: %d terms tested, %d with p.adjust<=0.1, top nominal: %s (p=%.4g)",
        label, name, nrow(df), n_sig, df$Description[1], df$pvalue[1])
    write_csv(utils::head(df, N_NOMINAL),
              file.path(outdir, sprintf("%s_%s_top%d_nominal.csv", label, name, N_NOMINAL)))
  }

  run_one(clusterProfiler::enrichGO, "GO_BP", OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
  run_one(clusterProfiler::enrichKEGG, "KEGG", organism = "hsa")
  run_one(ReactomePA::enrichPathway, "Reactome", organism = "human", readable = TRUE)
}


## ---- Main loop -------------------------------------------------------------------

for (run_name in RUNS) {
  outdir <- file.path(STEP6_DIR, run_name)
  cg_path <- file.path(outdir, "consensus_genes.csv")

  add(""); add("=====================================================================")
  add("RUN: %s", run_name)
  add("=====================================================================")

  if (!file.exists(cg_path)) { add("  !! consensus_genes.csv not found — skipping"); next }
  cg <- readr::read_csv(cg_path, show_col_types = FALSE)

  all_genes    <- unique(cg$gene)
  strict_genes <- cg |> dplyr::filter(tier != "extended (DIABLO only)") |> dplyr::pull(gene) |> unique()

  add("Full set: %d genes | Strict (core+moderate) set: %d genes", length(all_genes), length(strict_genes))

  add(""); add("-- FULL set --")
  run_hallmark(all_genes, outdir, paste0(run_name, "_FULL"))
  run_enrichment_with_nominal(all_genes, outdir, paste0(run_name, "_FULL"))

  if (length(strict_genes) >= 3) {
    add(""); add("-- STRICT set --")
    run_hallmark(strict_genes, outdir, paste0(run_name, "_STRICT"))
    run_enrichment_with_nominal(strict_genes, outdir, paste0(run_name, "_STRICT"))
  } else {
    add("  !! strict set too small (%d) for separate enrichment", length(strict_genes))
  }
}

REPORT_OUT <- file.path(STEP6_DIR, "Step6_1b_hallmark_nominal_report.txt")
writeLines(report, REPORT_OUT)
cat(paste(report, collapse = "\n"))
cat("\n\nFull report saved to:\n", REPORT_OUT, "\n")