# STEP 9.0 LOCAL AUDIT (no internet)
project_root <- "C:/Users/USER/Documents/MyThesis_MultiOmics"
out_dir <- file.path(project_root, "09_FinalResults", "Step9_Clinical")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
P <- function(...) file.path(project_root, ...)
short_id <- function(x) substr(as.character(x), 1, 12)
is_missing <- function(x) {
  x <- tolower(trimws(as.character(x)))
  is.na(x) | x %in% c("", "na", "nan", "#n/a", "not reported", "unknown",
    "[not available]", "[not evaluated]", "[unknown]", "[discrepancy]", "not applicable")
}
domains <- c(age = "^age|age_at", stage = "stage", grade = "grade|gleason",
             nodal = "pathologic_n|lymph|node", subtype = "subtype|pam50|molecular|histolog",
             survival = "days_to|vital|follow|death|recurr|progress|dfs|os_|_os")

files <- list(
  BRCA = c(P("01_RawData","Breast","clinical_brca.csv"),
           P("01_RawData","Breast","brca_rna_clinical_full.rds"),
           P("02_ProcessedData","Breast","brca_meth_clinical_full.rds")),
  PRAD = c(P("01_RawData","Prostate","clinical_prad.csv")))

load_tbl <- function(f) {
  if (!file.exists(f)) { cat("  MISSING:", f, "\n"); return(NULL) }
  x <- tryCatch(if (grepl("\\.rds$", f)) readRDS(f) else read.csv(f, check.names = FALSE, stringsAsFactors = FALSE),
                error = function(e) { cat("  read failed:", conditionMessage(e), "\n"); NULL })
  if (is.null(x)) return(NULL)
  if (inherits(x, "SummarizedExperiment")) x <- as.data.frame(SummarizedExperiment::colData(x))
  if (!is.data.frame(x)) { cat("  not a table; class:", class(x)[1], "\n"); return(NULL) }
  as.data.frame(x[, vapply(x, is.atomic, logical(1)), drop = FALSE])
}

for (p in names(files)) {
  d7 <- readRDS(P("07_FinalResults","Step7_FeatureSelection_ML", paste0(p, "_targeted_dataset.rds")))
  pid <- unique(short_id(d7$patient))
  cat("\n####################", p, "- Step 7 patients:", length(pid), "####################\n")
  for (f in files[[p]]) {
    cat("\n=== ", sub(project_root, "", f, fixed = TRUE), " ===\n", sep = "")
    x <- load_tbl(f); if (is.null(x)) next
    cat("rows:", nrow(x), "| cols:", ncol(x), "\n")
    # ID column with best match to Step 7 patients
    m <- vapply(names(x), function(cl) sum(pid %in% short_id(x[[cl]])), numeric(1))
    best <- names(which.max(m))
    cat("Best ID column:", best, "-> matches", max(m), "of", length(pid), "Step 7 patients\n")
    cat("Example IDs:", paste(head(x[[best]], 2), collapse = ", "), "\n")
    for (dm in names(domains)) {
      cols <- grep(domains[dm], names(x), value = TRUE, ignore.case = TRUE)
      if (length(cols) == 0) next
      cat("\n[", dm, "]\n", sep = "")
      for (cl in cols) cat(sprintf("  %-42s levels=%-4d missing=%5.1f%%\n", cl,
          length(unique(x[[cl]][!is_missing(x[[cl]])])), 100 * mean(is_missing(x[[cl]]))))
    }
    if (all(c("vital_status","days_to_death","days_to_last_follow_up") %in% names(x))) {
      dead <- as.integer(x$vital_status %in% c("Dead","DECEASED","dead"))
      t <- as.numeric(ifelse(dead == 1, x$days_to_death, x$days_to_last_follow_up))
      ok <- !is.na(t) & t > 0
      cat("\nOVERALL SURVIVAL: usable =", sum(ok), "| deaths =", sum(dead[ok]),
          "| median follow-up (yrs, crude) =", round(median(t[ok]) / 365.25, 2),
          "| max covariates @10 events each =", floor(sum(dead[ok]) / 10), "\n")
    }
  }
}

cat("\n#### FULL-COHORT RNA FILES ####\n")
for (f in c(P("02_ProcessedData","Breast","brca_vst.rds"), P("02_ProcessedData","Prostate","prad_rna_normalized.csv"),
            P("02_ProcessedData","Breast","brca_rna_normalized.csv"))) {
  cat("\n", sub(project_root, "", f, fixed = TRUE), "\n", sep = "")
  if (!file.exists(f)) { cat("  MISSING\n"); next }
  if (grepl("\\.rds$", f)) { x <- readRDS(f); cat("  class:", class(x)[1], "| dim:", paste(dim(x), collapse = " x "), "\n")
    if (!is.null(colnames(x))) cat("  first colnames:", paste(head(colnames(x), 3), collapse = ", "), "\n")
    if (!is.null(rownames(x))) cat("  first rownames:", paste(head(rownames(x), 3), collapse = ", "), "\n")
  } else { h <- names(read.csv(f, nrows = 1, check.names = FALSE))
    cat("  columns:", length(h), "| first:", paste(head(h, 4), collapse = ", "), "\n") }
}
cat("\nAudit done.\n")
