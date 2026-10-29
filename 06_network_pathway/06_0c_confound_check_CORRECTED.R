## =============================================================================
## STEP 6.0c — MOFA FACTOR1 CONFOUND CHECK (CORRECTED)
##
## The previous version (6.0b) used group ~ Factor1 + technical_metric via
## glm(binomial), which is invalid here: Factor1 separates the groups almost
## perfectly (Kruskal-Wallis p ~1e-6), causing logistic-regression quasi-
## separation — coefficients diverge and Wald p-values become meaningless
## (inflated toward 1), producing a false "not associated" result.
##
## Fix: residualize Factor1 against the technical metric via linear
## regression, then test the RESIDUAL against group with Kruskal-Wallis
## (same non-parametric test used to originally select Factor1 — no
## separation problem). This tells us whether group-discriminating signal
## survives once the technical component is removed from Factor1.
## =============================================================================

library(dplyr)
library(readr)
library(tibble)

ROOT      <- "C:/Users/USER/Documents/MyThesis_MultiOmics"
BASE      <- file.path(ROOT, "05_FinalResults", "Step5.0_Integration")
STEP6_DIR <- file.path(ROOT, "06_FinalResults", "Step6_Network")

RUNS <- list(
  BRCA_3omics   = list(dir = file.path(BASE, "BRCA", "3omics_16pairs"),   prep = "A3"),
  BRCA_meth_rna = list(dir = file.path(BASE, "BRCA", "meth_rna_51pairs"), prep = "A2"),
  PRAD_meth_rna = list(dir = file.path(BASE, "PRAD", "meth_rna_17pairs"), prep = "A")
)

report <- c()
add <- function(...) report <<- c(report, sprintf(...))

per_sample_technical_metric <- function(mat, sample_ids, block_label) {
  if (is.null(mat)) return(NULL)
  n <- length(sample_ids)
  if (ncol(mat) == n) {
    metric <- colMeans(mat, na.rm = TRUE); ids_used <- colnames(mat)
  } else if (nrow(mat) == n) {
    metric <- rowMeans(mat, na.rm = TRUE); ids_used <- rownames(mat)
  } else {
    add("  !! block '%s' dims don't match sample count — skipped", block_label)
    return(NULL)
  }
  tibble(sample_id = ids_used, metric = as.numeric(metric), block = block_label)
}

for (run_name in names(RUNS)) {
  cfg <- RUNS[[run_name]]; dir_path <- cfg$dir
  outdir <- file.path(STEP6_DIR, run_name)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  add(""); add("=====================================================================")
  add("RUN: %s", run_name)
  add("=====================================================================")

  e <- new.env(); load(file.path(dir_path, "..", "prep_objects.RData"), envir = e)
  prep <- get(cfg$prep, envir = e)
  X <- prep$X; group_vec <- prep$group; sample_ids <- prep$sid

  tech_metrics <- list()
  for (blk in names(X)) {
    tm <- per_sample_technical_metric(X[[blk]], sample_ids, blk)
    if (!is.null(tm)) tech_metrics[[blk]] <- tm
  }

  hdf5_path <- file.path(dir_path, "MOFA_model.hdf5")
  model <- tryCatch(MOFA2::load_model(hdf5_path), error = function(e) NULL)
  if (is.null(model)) { add("  !! failed to load MOFA model"); next }
  factors_df <- as.data.frame(MOFA2::get_factors(model)[[1]]) |> rownames_to_column("sample_id")
  meta <- tibble(sample_id = sample_ids, group = group_vec)
  base_df <- factors_df |> inner_join(meta, by = "sample_id")

  ## Find the block whose technical metric correlates most strongly with Factor1
  best_rho <- 0; best_block <- NA_character_
  for (blk in names(tech_metrics)) {
    merged <- base_df |> inner_join(tech_metrics[[blk]], by = "sample_id")
    rho <- suppressWarnings(cor(merged$Factor1, merged$metric, method = "spearman"))
    if (!is.na(rho) && abs(rho) > abs(best_rho)) { best_rho <- rho; best_block <- blk }
  }
  add("Strongest Factor1 <-> technical-metric correlation: block '%s', rho = %.3f", best_block, best_rho)

  merged <- base_df |> inner_join(tech_metrics[[best_block]], by = "sample_id")

  ## Original (raw) Factor1 vs group
  p_raw <- kruskal.test(merged$Factor1 ~ as.factor(merged$group))$p.value

  ## Residualize Factor1 against the technical metric, then test residual vs group
  resid_fit <- lm(Factor1 ~ metric, data = merged)
  merged$Factor1_resid <- resid(resid_fit)
  p_resid <- kruskal.test(merged$Factor1_resid ~ as.factor(merged$group))$p.value

  r2_explained <- summary(resid_fit)$r.squared

  add("")
  add("  Factor1 (raw)      vs group : Kruskal-Wallis p = %s", signif(p_raw, 4))
  add("  Factor1 (residual) vs group : Kruskal-Wallis p = %s", signif(p_resid, 4))
  add("  Variance in Factor1 explained by technical metric: R^2 = %.3f", r2_explained)

  verdict <- if (p_resid < 0.05) {
    sprintf(paste("Factor1 remains significantly associated with the group label (p = %s) after removing",
                   "the technical-metric component (which explained %.1f%% of Factor1's variance).",
                   "Group-discriminating signal survives independent of the normalization confound,",
                   "though %.1f%% of Factor1 is attributable to the technical effect and should be disclosed."),
                   signif(p_resid, 4), r2_explained * 100, r2_explained * 100)
  } else {
    sprintf(paste("Factor1's association with group does not survive after removing the technical-metric",
                   "component (residual p = %s vs raw p = %s), and the technical metric explains %.1f%% of",
                   "Factor1's variance. This suggests Factor1's group-discriminating power may be substantially",
                   "confounded by normalization/library-size effects. Recommend treating MOFA support from",
                   "Factor1 with caution, or excluding it from the consensus MOFA-support criterion for this run."),
                   signif(p_resid, 4), signif(p_raw, 4), r2_explained * 100)
  }
  add(""); add("  VERDICT: %s", verdict)

  write_csv(
    tibble(run = run_name, technical_block = best_block, factor1_rho_with_technical = best_rho,
           r2_explained_by_technical = r2_explained, p_raw = p_raw, p_residual = p_resid, verdict = verdict),
    file.path(outdir, "factor1_confound_verdict_CORRECTED.csv")
  )
}

REPORT_OUT <- file.path(STEP6_DIR, "Step6_0c_confound_check_CORRECTED_report.txt")
writeLines(report, REPORT_OUT)
cat(paste(report, collapse = "\n"))
cat("\n\nFull report saved to:\n", REPORT_OUT, "\n")
