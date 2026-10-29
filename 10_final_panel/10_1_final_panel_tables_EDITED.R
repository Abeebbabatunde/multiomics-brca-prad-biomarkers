# EDITED COPY: only the hard-coded protein sample size text differs (n=76 -> n=30, RPPA). Everything else is unchanged.
# STEP 10.1-10.3 FINAL PANEL TABLES - PART 1/2
project_root <- "C:/Users/USER/Documents/MyThesis_MultiOmics"
P <- function(...) file.path(project_root, ...)
out_dir <- P("10_FinalResults", "Step10_FinalPanel")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
s9 <- P("09_FinalResults", "Step9_Clinical")
s8 <- P("08_FinalResults", "Step8_External_Validation")

rd_csv <- function(f) {
  if (!file.exists(f)) stop("Missing file: ", f, call. = FALSE)
  read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
}

step7 <- readRDS(P("07_FinalResults", "Step7_FeatureSelection_ML",
                   "Step7_feature_selection_results.rds"))
panel <- list(BRCA = step7$BRCA$consensus, PRAD = step7$PRAD$consensus)
protein_feats <- list(BRCA = "VAV1", PRAD = character(0))
full <- list(
  BRCA = readRDS(P("07_FinalResults", "Step7_FeatureSelection_ML", "BRCA_targeted_dataset.rds")),
  PRAD = readRDS(P("07_FinalResults", "Step7_FeatureSelection_ML", "PRAD_targeted_dataset.rds")))

ext  <- rd_csv(file.path(s8, "Step8_FINAL_external_summary.csv"))
clin <- rd_csv(file.path(s9, "Step9_1_clinical_association_results.csv"))
adj  <- rd_csv(file.path(s9, "Step9_2_cox_multivariable.csv"))
km   <- rd_csv(file.path(s9, "Step9_2_KM_logrank.csv"))

omics_of <- function(cancer, f) {
  ifelse(grepl("^cg", f), "methylation",
         ifelse(f %in% protein_feats[[cancer]], "protein", "RNA"))
}

# descriptive direction from the paired Step 7 data (units as stored there)
paired_direction <- function(d, feats) {
  d <- as.data.frame(d)
  do.call(rbind, lapply(feats, function(f) {
    out <- data.frame(feature = f, n_pairs = 0L, median_diff = NA_real_,
                      p_paired = NA_real_, direction = NA_character_,
                      stringsAsFactors = FALSE)
    if (!f %in% names(d)) return(out)
    a <- d[d$group == "Normal", c("patient", f)]
    b <- d[d$group == "Early",  c("patient", f)]
    a <- a[!duplicated(a$patient), ]; b <- b[!duplicated(b$patient), ]
    names(a)[2] <- "normal"; names(b)[2] <- "early"
    m <- merge(a, b, by = "patient")
    m$normal <- as.numeric(m$normal); m$early <- as.numeric(m$early)
    m <- m[complete.cases(m), ]
    out$n_pairs <- nrow(m)
    if (nrow(m) >= 5) {
      dd <- m$early - m$normal
      out$median_diff <- median(dd)
      out$p_paired <- wilcox.test(m$early, m$normal, paired = TRUE, exact = FALSE)$p.value
      out$direction <- ifelse(median(dd) > 0, "higher in tumor", "lower in tumor")
    }
    out
  }))
}

ext_phrase <- function(cancer, f, om) {
  if (om == "protein") return("not evaluated (no independent external RPPA cohort)")
  r <- ext[ext$cohort == cancer & ext$omics == om, , drop = FALSE]
  if (!nrow(r) || is.na(r$auc[1])) return("external layer model not available")
  used <- strsplit(as.character(r$features_used[1]), ";")[[1]]
  base <- sprintf("%s layer model AUC %.3f (%s, n=%d)", om, r$auc[1], r$geo[1], r$n[1])
  if (f %in% used) paste0(base, "; feature included")
  else paste0(base, "; feature NOT evaluable (absent from external platform)")
}

sig_clin <- function(cancer, f) {
  r <- clin[clin$cancer == cancer & clin$feature == f & !is.na(clin$fdr) & clin$fdr < 0.05, ,
            drop = FALSE]
  if (!nrow(r)) return(NA_character_)
  paste(sprintf("%s (%s=%.2f, FDR=%.2g)", r$variable, r$effect_measure, r$effect, r$fdr),
        collapse = "; ")
}

prog_info <- function(cancer, f) {
  if (cancer == "PRAD")
    return(list(txt = "not evaluable (2 deaths among 188 patients)",
                status = "prognostic evidence unavailable"))
  a <- adj[adj$feature == f, , drop = FALSE]
  if (!nrow(a) || is.na(a$fdr[1]))
    return(list(txt = "not evaluated (protein, n=30, too few events)",
                status = "prognostic evidence unavailable"))
  k <- km[km$feature == f, , drop = FALSE]
  kfdr <- if (nrow(k) && !is.na(k$fdr[1])) sprintf("%.2f", k$fdr[1]) else "NA"
  list(txt = sprintf("adjusted HR per SD %.2f (%.2f-%.2f), FDR=%.2f; KM log-rank FDR=%s",
                     a$HR_per_SD[1], a$ci_low[1], a$ci_high[1], a$fdr[1], kfdr),
       status = if (a$fdr[1] < 0.05) "prognostic association supported"
                else "prognostic association not established")
}

build_panel <- function(cancer) {
  feats <- panel[[cancer]]
  dir <- paired_direction(full[[cancer]], feats)
  do.call(rbind, lapply(feats, function(f) {
    om <- omics_of(cancer, f)
    ex <- ext_phrase(cancer, f, om)
    cs <- sig_clin(cancer, f)
    pg <- prog_info(cancer, f)
    dr <- dir[dir$feature == f, ]
    clin_status <- if (is.na(cs)) "no clinical association established" else "clinically associated"
    needs <- grepl("NOT evaluable|not evaluated", ex)
    data.frame(
      cancer = cancer, feature = f, omics = om,
      tumor_vs_normal = dr$direction, n_pairs = dr$n_pairs,
      median_paired_diff = signif(dr$median_diff, 3), paired_p = signif(dr$p_paired, 3),
      external_classification_evidence = ex,
      clinical_evidence_step9_1 = ifelse(is.na(cs), "none at FDR<0.05", cs),
      prognostic_evidence_step9_2 = pg$txt,
      clinical_status = clin_status, prognostic_status = pg$status,
      needs_further_validation = ifelse(needs, "external classification validation", "no"),
      final_interpretation = paste0("Retained (locked Step 7); ", clin_status, "; ", pg$status,
        ifelse(needs, "; requires further validation (external classification)", "")),
      stringsAsFactors = FALSE)
  }))
}

# PART 2/2

panels_out <- list()
for (cn in c("BRCA", "PRAD")) {
  tb <- build_panel(cn); panels_out[[cn]] <- tb
  fn <- paste0("Step10_", ifelse(cn == "BRCA", "1_BRCA", "2_PRAD"), "_final_panel.csv")
  write.csv(tb, file.path(out_dir, fn), row.names = FALSE)
  cat("\n=====", cn, "FINAL PANEL:", nrow(tb), "features =====\n")
  print(tb[, c("feature", "omics", "tumor_vs_normal", "clinical_status",
               "prognostic_status", "needs_further_validation")], row.names = FALSE)
}

# ---- 10.3 panel-level comparison
internal <- readRDS(P("08_FinalResults", "Step8_Validation", "Step8_internal_validation_summary.rds"))
get_int <- function(cn, what) {
  v <- tryCatch({
    s <- internal[[cn]]
    if (what == "cv") { cv <- s$cv_performance; cv$mean[cv$metric == "AUC"] } else s$holdout_test$auc
  }, error = function(e) NA_real_)
  if (length(v) == 0) NA_real_ else as.numeric(v[1])
}
auc_of <- function(cn, om) {
  v <- ext$auc[ext$cohort == cn & ext$omics == om]
  if (length(v)) v[1] else NA_real_
}
cpanel <- tryCatch(rd_csv(file.path(s9, "Step9_2_RNA_panel_vs_clinical.csv")), error = function(e) NULL)

perf <- do.call(rbind, lapply(c("BRCA", "PRAD"), function(cn) {
  f <- panel[[cn]]; om <- omics_of(cn, f); tb <- panels_out[[cn]]
  data.frame(cancer = cn, n_features = length(f),
    n_methylation = sum(om == "methylation"), n_RNA = sum(om == "RNA"),
    n_protein = sum(om == "protein"),
    internal_cv_auc = get_int(cn, "cv"), internal_holdout_auc = get_int(cn, "holdout"),
    external_methylation_auc = auc_of(cn, "methylation"),
    external_RNA_auc = auc_of(cn, "RNA"),
    features_clinically_associated = sum(tb$clinical_status == "clinically associated"),
    features_prognostic_supported = sum(tb$prognostic_status == "prognostic association supported"),
    features_needing_validation = sum(tb$needs_further_validation != "no"),
    rna_panel_C_clinical_only = if (cn == "BRCA" && !is.null(cpanel)) cpanel$C_clinical_only[1] else NA_real_,
    rna_panel_C_with_rna = if (cn == "BRCA" && !is.null(cpanel)) cpanel$C_clinical_plus_RNA[1] else NA_real_,
    rna_panel_LR_p = if (cn == "BRCA" && !is.null(cpanel)) cpanel$LR_p_value[1] else NA_real_,
    stringsAsFactors = FALSE)
}))
write.csv(perf, file.path(out_dir, "Step10_3_panel_comparison.csv"), row.names = FALSE)

shared <- intersect(panel$BRCA, panel$PRAD)
non_meth <- function(cn) panel[[cn]][omics_of(cn, panel[[cn]]) != "methylation"]
shared_genes <- intersect(non_meth("BRCA"), non_meth("PRAD"))
overlap <- data.frame(
  level = c("exact feature (probe or gene)", "RNA/protein gene symbol"),
  n_shared = c(length(shared), length(shared_genes)),
  shared = c(paste(shared, collapse = ";"), paste(shared_genes, collapse = ";")))
write.csv(overlap, file.path(out_dir, "Step10_3_feature_overlap.csv"), row.names = FALSE)

cat("\n=== 10.3 PANEL COMPARISON ===\n"); print(t(perf))
cat("\n=== FEATURE OVERLAP ===\n"); print(overlap, row.names = FALSE)
cat("\nSaved in:", out_dir, "\n")
