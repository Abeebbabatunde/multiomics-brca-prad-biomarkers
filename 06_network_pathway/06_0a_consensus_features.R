## =============================================================================
## STEP 6.0a — CONSENSUS FEATURE BUILDER
## Builds the tiered DIABLO + MOFA + univariate-DE consensus gene set for each
## of the three real integration runs, using the actual files/objects found:
##   BRCA/3omics_16pairs, BRCA/meth_rna_51pairs, PRAD/meth_rna_17pairs
## Output: one consensus table per run, ready to feed into the Step 6 network
## /pathway script (STEP6_DIR/<run>/consensus_genes.csv).
## =============================================================================

library(dplyr)
library(readr)
library(stringr)
library(tibble)

ROOT      <- "C:/Users/USER/Documents/MyThesis_MultiOmics"
BASE      <- file.path(ROOT, "05_FinalResults", "Step5.0_Integration")
STEP6_DIR <- file.path(ROOT, "06_FinalResults", "Step6_Network")

RUNS <- list(
  BRCA_3omics   = list(dir = file.path(BASE, "BRCA", "3omics_16pairs"),   prep = "A3"),
  BRCA_meth_rna = list(dir = file.path(BASE, "BRCA", "meth_rna_51pairs"), prep = "A2"),
  PRAD_meth_rna = list(dir = file.path(BASE, "PRAD", "meth_rna_17pairs"), prep = "A")
)

## PRAD fallback DE files (since signature_with_DE.csv doesn't exist for PRAD)
PRAD_DMP_FILE <- file.path(ROOT, "05_FinalResults", "..", "03_Results", "Epigenomics", "prad_DMPs_full_results.csv")
PRAD_DEG_FILE <- file.path(ROOT, "05_FinalResults", "..", "03_Results", "Transcriptomics", "prad_DEGs_full_results.csv")
# NOTE: adjust the above two paths if your actual project root differs from ROOT/03_Results

report <- c()
add <- function(...) report <<- c(report, sprintf(...))


## ---- helper: CpG -> gene symbol (450K) ---------------------------------------
cpg_to_symbol <- function(cpg_ids) {
  suppressMessages({
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    library(minfi)
  })
  ann <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
  hit <- ann[cpg_ids, c("UCSC_RefGene_Name"), drop = FALSE]
  genes <- sapply(strsplit(hit$UCSC_RefGene_Name, ";"), function(x) if (length(x) == 0) NA else x[1])
  tibble(feature = cpg_ids, gene = genes)
}


## ---- helper: identify the label-associated MOFA factor -----------------------
identify_label_factor <- function(hdf5_path, group_vec, sample_ids) {
  out <- list(factor = NA_character_, pvals = NULL, matched_n = 0)
  if (!file.exists(hdf5_path)) { add("  !! MOFA hdf5 not found: %s", hdf5_path); return(out) }
  ok <- requireNamespace("MOFA2", quietly = TRUE)
  if (!ok) { add("  !! MOFA2 package not available — skipping factor-label test"); return(out) }

  model <- tryCatch(MOFA2::load_model(hdf5_path), error = function(e) {
    add("  !! load_model failed: %s", conditionMessage(e)); NULL
  })
  if (is.null(model)) return(out)

  factors_list <- MOFA2::get_factors(model)
  factors_df <- as.data.frame(factors_list[[1]]) |> rownames_to_column("mofa_sample")

  meta <- tibble(sample_id = sample_ids, group = group_vec)
  merged <- factors_df |> inner_join(meta, by = c("mofa_sample" = "sample_id"))
  out$matched_n <- nrow(merged)
  add("  MOFA samples: %d | group vec length: %d | matched by name: %d",
      nrow(factors_df), length(sample_ids), nrow(merged))

  if (nrow(merged) < 4) {
    add("  !! too few matched samples to test factor<->label association (name mismatch likely)")
    return(out)
  }

  merged$group <- as.factor(merged$group)
  fac_cols <- setdiff(names(factors_df), "mofa_sample")
  pvals <- sapply(fac_cols, function(fc) {
    tryCatch(kruskal.test(merged[[fc]] ~ merged$group)$p.value, error = function(e) NA)
  })
  add("  Factor <-> group association (Kruskal-Wallis p):")
  for (fc in fac_cols) add("    %s : p = %s", fc, ifelse(is.na(pvals[fc]), "NA", signif(pvals[fc], 4)))

  best <- names(which.min(pvals))
  out$factor <- best
  out$pvals <- pvals
  add("  --> selected label-associated factor: %s (p = %s)", best, signif(min(pvals, na.rm = TRUE), 4))
  out
}


## ---- main loop per run --------------------------------------------------------

all_consensus <- list()

for (run_name in names(RUNS)) {
  cfg <- RUNS[[run_name]]
  dir_path <- cfg$dir
  outdir <- file.path(STEP6_DIR, run_name)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  add("")
  add("=====================================================================")
  add("RUN: %s", run_name)
  add("=====================================================================")

  ## 1. DIABLO selected features (all blocks)
  diablo <- read_csv(file.path(dir_path, "DIABLO_selected_features.csv"), show_col_types = FALSE)
  add("DIABLO selected features by block:")
  print(table(diablo$block))
  add("  total selected: %d", nrow(diablo))

  ## 2. Load prep_objects.RData for group labels + sample IDs
  e <- new.env()
  load(file.path(dir_path, "..", "prep_objects.RData"), envir = e)  # one level up (BRCA/ or PRAD/)
  prep <- get(cfg$prep, envir = e)
  group_vec <- prep$group
  sample_ids <- prep$sid
  add("Group label distribution:")
  print(table(group_vec))

  ## 3. Identify label-associated MOFA factor
  add("")
  add("MOFA factor identification:")
  hdf5_path <- file.path(dir_path, "MOFA_model.hdf5")
  factor_info <- identify_label_factor(hdf5_path, group_vec, sample_ids)

  ## 4. Load MOFA weights for the identified factor (all views), threshold top |value|
  weight_files <- list.files(dir_path, pattern = "^MOFA_weights_.*\\.csv$", full.names = TRUE)
  mofa_top <- tibble()
  if (!is.na(factor_info$factor)) {
    for (wf in weight_files) {
      wdf <- read_csv(wf, show_col_types = FALSE) |> filter(factor == factor_info$factor)
      # top 1% by absolute weight, per view
      cutoff <- quantile(abs(wdf$value), 0.99, na.rm = TRUE)
      top <- wdf |> filter(abs(value) >= cutoff)
      mofa_top <- bind_rows(mofa_top, top)
    }
    add("MOFA top-loading features on %s (top 1%% by |weight| per view): %d total",
        factor_info$factor, nrow(mofa_top))
  } else {
    add("  !! no label-associated factor identified — MOFA support will be skipped for this run")
  }

  ## 5. DE overlap (signature_with_DE.csv, or built on the fly for PRAD)
  sig_path <- file.path(dir_path, "signature_with_DE.csv")
  if (file.exists(sig_path)) {
    sig_de <- read_csv(sig_path, show_col_types = FALSE)
  } else {
    add("  signature_with_DE.csv missing — attempting to build from PRAD DMP/DEG result files")
    sig_de <- diablo  # start from DIABLO panel
    sig_de$de_effect <- NA_real_; sig_de$de_padj <- NA_real_; sig_de$sig_in_DE <- FALSE

    if (file.exists(PRAD_DMP_FILE)) {
      dmp <- read_csv(PRAD_DMP_FILE, show_col_types = FALSE)
      id_col <- names(dmp)[1]  # assume first column is the CpG/probe ID
      padj_col <- names(dmp)[str_detect(tolower(names(dmp)), "padj|adj.*p|fdr")][1]
      meth_rows <- sig_de$block == "meth"
      match_idx <- match(sig_de$feature[meth_rows], dmp[[id_col]])
      sig_de$de_padj[meth_rows] <- dmp[[padj_col]][match_idx]
      sig_de$sig_in_DE[meth_rows] <- !is.na(sig_de$de_padj[meth_rows]) & sig_de$de_padj[meth_rows] < 0.05
    } else {
      add("  !! PRAD_DMP_FILE not found at %s — edit path at top of script", PRAD_DMP_FILE)
    }
    if (file.exists(PRAD_DEG_FILE)) {
      deg <- read_csv(PRAD_DEG_FILE, show_col_types = FALSE)
      id_col <- names(deg)[1]
      padj_col <- names(deg)[str_detect(tolower(names(deg)), "padj|adj.*p|fdr")][1]
      rna_rows <- sig_de$block == "rna"
      match_idx <- match(sig_de$feature[rna_rows], deg[[id_col]])
      sig_de$de_padj[rna_rows] <- deg[[padj_col]][match_idx]
      sig_de$sig_in_DE[rna_rows] <- !is.na(sig_de$de_padj[rna_rows]) & sig_de$de_padj[rna_rows] < 0.05
    } else {
      add("  !! PRAD_DEG_FILE not found at %s — edit path at top of script", PRAD_DEG_FILE)
    }
  }
  add("DE-overlap support: %d / %d DIABLO features sig_in_DE == TRUE",
      sum(sig_de$sig_in_DE, na.rm = TRUE), nrow(sig_de))

  ## 6. Build tiered consensus at the FEATURE level
  sig_de <- sig_de |>
    mutate(mofa_support = feature %in% mofa_top$feature)

  consensus <- sig_de |>
    mutate(
      n_support = as.integer(sig_in_DE) + as.integer(mofa_support),  # DIABLO always counts (it's the base panel)
      tier = case_when(
        n_support >= 2 ~ "core (DIABLO + DE + MOFA)",
        n_support == 1 ~ "moderate (DIABLO + 1 corroboration)",
        TRUE           ~ "extended (DIABLO only)"
      )
    )

  ## 7. Map to gene symbols (CpG -> gene; rna/prot already symbols) and collapse to gene level
  cpg_feats <- consensus |> filter(block == "meth") |> pull(feature)
  cpg_map <- if (length(cpg_feats) > 0) cpg_to_symbol(cpg_feats) else tibble(feature = character(), gene = character())

  consensus_genes <- consensus |>
    left_join(cpg_map, by = "feature") |>
    mutate(gene = ifelse(block == "meth", gene, feature)) |>
    filter(!is.na(gene))

  write_csv(consensus, file.path(outdir, "consensus_features_raw.csv"))
  write_csv(consensus_genes, file.path(outdir, "consensus_genes.csv"))

  add("Tier breakdown (feature level):")
  print(table(consensus$tier))
  add("Unique consensus genes written: %d -> %s", length(unique(consensus_genes$gene)),
      file.path(outdir, "consensus_genes.csv"))

  all_consensus[[run_name]] <- consensus_genes
}


## ---- Write report --------------------------------------------------------------

REPORT_OUT <- file.path(STEP6_DIR, "Step6_0a_consensus_report.txt")
writeLines(report, REPORT_OUT)
cat(paste(report, collapse = "\n"))
cat("\n\nFull report saved to:\n", REPORT_OUT, "\n")
cat("\nPer-run consensus_genes.csv files written under:\n", STEP6_DIR, "\n")
