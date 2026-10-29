# STEP 9.1 CLINICAL ASSOCIATION - PART 1/3
# EDITED COPY: only get_protein() differs from the run version (reads RPPA instead of the iTRAQ file). Everything else is unchanged.
project_root <- "C:/Users/USER/Documents/MyThesis_MultiOmics"
P <- function(...) file.path(project_root, ...)
out_dir <- P("09_FinalResults", "Step9_Clinical")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(42)
min_group_n <- 5      # smallest group allowed in a comparison
min_cont_n  <- 20     # smallest n for a correlation
n_boot      <- 1000   # bootstrap resamples for epsilon-squared CIs
manual_ensembl <- character(0)  # e.g. c(CHD9NB = "ENSG...") if a symbol fails to map

for (pk in c("data.table", "org.Hs.eg.db", "AnnotationDbi",
             "SummarizedExperiment", "DESeq2")) {
  if (!requireNamespace(pk, quietly = TRUE))
    stop("Install this package first: ", pk, call. = FALSE)
}

step7 <- readRDS(P("07_FinalResults", "Step7_FeatureSelection_ML",
                   "Step7_feature_selection_results.rds"))
panel <- list(BRCA = step7$BRCA$consensus, PRAD = step7$PRAD$consensus)
protein_feats <- list(BRCA = "VAV1", PRAD = character(0))

omics_of <- function(cancer, feats) {
  ifelse(grepl("^cg", feats), "methylation",
         ifelse(feats %in% protein_feats[[cancer]], "protein", "RNA"))
}

short_id <- function(x) substr(as.character(x), 1, 12)

is_missing <- function(x) {
  x <- tolower(trimws(as.character(x)))
  is.na(x) | x %in% c("", "na", "nan", "not reported", "unknown",
                      "[not available]", "[unknown]", "[not evaluated]")
}

# one primary-tumor sample (type 01) per patient; first barcode in sort order
one_tumor_per_patient <- function(samples) {
  s <- sort(samples[substr(samples, 14, 15) == "01"])
  s[!duplicated(short_id(s))]
}

# read selected rows of a big CSV without loading the whole file
pull_rows_csv <- function(f, wanted, strip_version = FALSE) {
  ids <- data.table::fread(f, select = 1, header = TRUE)[[1]]
  key <- if (strip_version) sub("\\..*$", "", ids) else ids
  idx <- match(wanted, key)
  hdr <- names(data.table::fread(f, nrows = 0))
  out <- matrix(NA_real_, length(wanted), length(hdr) - 1,
                dimnames = list(as.character(wanted), hdr[-1]))
  for (i in which(!is.na(idx))) {
    r <- as.data.frame(data.table::fread(f, skip = idx[i], nrows = 1, header = FALSE))
    out[i, ] <- as.numeric(unlist(r[1, -1]))
  }
  out
}

map_to_ensembl <- function(symbols, present_ids) {
  res <- setNames(rep(NA_character_, length(symbols)), symbols)
  for (s in symbols) {
    if (s %in% names(manual_ensembl)) { res[s] <- manual_ensembl[[s]]; next }
    for (kt in c("SYMBOL", "ALIAS")) {
      m <- tryCatch(suppressMessages(AnnotationDbi::select(
             org.Hs.eg.db::org.Hs.eg.db, keys = s, keytype = kt, columns = "ENSEMBL")),
           error = function(e) NULL)
      if (is.null(m)) next
      cand <- unique(stats::na.omit(m$ENSEMBL))
      cand <- cand[cand %in% present_ids]
      if (length(cand) > 0) {
        res[s] <- cand[1]
        if (kt == "ALIAS") cat("  NOTE:", s, "mapped via alias -> check", cand[1], "\n")
        break
      }
    }
    cat(sprintf("  %-8s -> %s\n", s, ifelse(is.na(res[s]), "NOT FOUND", res[s])))
  }
  res
}

# PART 2/3

build_clinical <- function(cancer) {
  if (cancer == "BRCA") {
    cb <- read.csv(P("01_RawData", "Breast", "clinical_brca.csv"),
                   check.names = FALSE, stringsAsFactors = FALSE)
    st <- cb$ajcc_pathologic_stage
    nn <- cb$ajcc_pathologic_n
    cl <- data.frame(
      patient = cb$submitter_id,
      age     = as.numeric(cb$age_at_index),
      stage   = ifelse(grepl("^Stage I[AB]?$", st), "Stage I",
                ifelse(grepl("^Stage II[AB]?$", st), "Stage II", NA)),
      nodal   = ifelse(grepl("^N0", nn), "N0",
                ifelse(grepl("^N[1-3]", nn), "N+", NA)),
      stringsAsFactors = FALSE)
    cl$stage <- factor(cl$stage, levels = c("Stage I", "Stage II"))
    cl$nodal <- factor(cl$nodal, levels = c("N0", "N+"))
    mc <- as.data.frame(readRDS(P("02_ProcessedData", "Breast", "brca_meth_clinical_full.rds")))
    pam <- data.frame(patient = short_id(mc$barcode),
                      type    = substr(mc$barcode, 14, 15),
                      pam50   = as.character(mc$paper_BRCA_Subtype_PAM50),
                      stringsAsFactors = FALSE)
    pam <- pam[pam$type == "01" & !is_missing(pam$pam50), ]
    pam <- pam[!duplicated(pam$patient), ]
    cl$subtype <- factor(pam$pam50[match(cl$patient, pam$patient)])
  } else {
    cp <- read.csv(P("01_RawData", "Prostate", "clinical_prad.csv"),
                   check.names = FALSE, stringsAsFactors = FALSE)
    gs <- suppressWarnings(as.numeric(cp$gleason_score))
    nn <- trimws(cp$ajcc_pathologic_n)
    cl <- data.frame(
      patient   = cp$submitter_id,
      age       = as.numeric(cp$age_at_index),
      gleason   = ifelse(is.na(gs), NA, ifelse(gs <= 6, "Gleason 6",
                  ifelse(gs == 7, "Gleason 7", "Gleason 8-9"))),
      highgrade = ifelse(is.na(gs), NA, ifelse(gs >= 8, "Gleason 8-9", "Gleason 6-7")),
      nodal     = ifelse(nn == "N0", "N0", ifelse(nn == "N1", "N1", NA)),
      stringsAsFactors = FALSE)
    cl$gleason   <- factor(cl$gleason, levels = c("Gleason 6", "Gleason 7", "Gleason 8-9"))
    cl$highgrade <- factor(cl$highgrade, levels = c("Gleason 6-7", "Gleason 8-9"))
    cl$nodal     <- factor(cl$nodal, levels = c("N0", "N1"))
  }
  cl
}

to_df <- function(vals) {
  df <- data.frame(patient = short_id(rownames(vals)), vals,
                   check.names = FALSE, stringsAsFactors = FALSE)
  rownames(df) <- NULL
  df[!duplicated(df$patient), ]
}

get_rna <- function(cancer, genes) {
  if (length(genes) == 0) return(NULL)
  if (cancer == "BRCA") {
    v <- readRDS(P("02_ProcessedData", "Breast", "brca_vst.rds"))
    mat <- SummarizedExperiment::assay(v)
    ids <- sub("\\..*$", "", rownames(mat))
    ens <- map_to_ensembl(genes, ids)
    keep <- names(ens)[!is.na(ens)]
    if (!length(keep)) return(NULL)
    tum <- one_tumor_per_patient(colnames(mat))
    vals <- t(mat[match(ens[keep], ids), tum, drop = FALSE])
  } else {
    f <- P("02_ProcessedData", "Prostate", "prad_rna_normalized.csv")
    ids <- sub("\\..*$", "", data.table::fread(f, select = 1)[[1]])
    ens <- map_to_ensembl(genes, ids)
    keep <- names(ens)[!is.na(ens)]
    if (!length(keep)) return(NULL)
    m <- pull_rows_csv(f, ens[keep], strip_version = TRUE)
    tum <- one_tumor_per_patient(colnames(m))
    vals <- t(m[, tum, drop = FALSE])
  }
  colnames(vals) <- keep
  to_df(vals)
}

get_meth <- function(cancer, probes) {
  if (!length(probes)) return(NULL)
  f <- if (cancer == "BRCA") P("02_ProcessedData", "Breast", "brca_meth_normalized.csv")
       else P("02_ProcessedData", "Prostate", "prad_meth_normalized.csv")
  m <- pull_rows_csv(f, probes)
  miss <- probes[rowSums(!is.na(m)) == 0]
  if (length(miss)) cat("  probes not found:", paste(miss, collapse = ", "), "\n")
  m <- m[rowSums(!is.na(m)) > 0, , drop = FALSE]
  if (nrow(m) == 0) return(NULL)
  tum <- one_tumor_per_patient(colnames(m))
  to_df(t(m[, tum, drop = FALSE]))
}

get_protein <- function(cancer, feats) {
  if (!length(feats)) return(NULL)
  # EDITED: protein values come from the RPPA layer used in the rest of the thesis
  # (long format: peptide_target, sample_barcode, value, sample_type), tumour samples only.
  pr <- data.table::fread(P("02_Normalized", "Breast", "brca_protein_normalized.csv"))
  pr <- pr[sample_type == "Primary Tumor" & peptide_target %in% feats]
  hit <- feats[feats %in% pr$peptide_target]
  if (!length(hit)) { cat("  protein feature not found\n"); return(NULL) }
  w <- data.table::dcast(pr, sample_barcode ~ peptide_target, value.var = "value")
  vals <- as.matrix(w[, hit, with = FALSE])
  rownames(vals) <- w$sample_barcode
  to_df(vals)
}

build_dataset <- function(cancer) {
  cat("\n######", cancer, "######\n")
  clin  <- build_clinical(cancer)
  feats <- panel[[cancer]]
  om    <- omics_of(cancer, feats)
  parts <- list(get_rna(cancer, feats[om == "RNA"]),
                get_meth(cancer, feats[om == "methylation"]),
                get_protein(cancer, feats[om == "protein"]))
  parts <- Filter(Negate(is.null), parts)
  fdf   <- Reduce(function(a, b) merge(a, b, by = "patient", all = TRUE), parts)
  dat   <- merge(clin, fdf, by = "patient", all.x = TRUE)
  fc    <- setdiff(names(fdf), "patient")
  cat("Clinical patients:", nrow(clin), "| features available:", length(fc),
      "of", length(feats), "\n")
  list(data = dat, clinical = clin,
       features = data.frame(feature = fc, omics = omics_of(cancer, fc),
                             stringsAsFactors = FALSE))
}

# PART 3/3 (revised: Cliff's delta for two-group comparisons)
set.seed(42)

# positive = higher in the second factor level
cliff_delta <- function(x, g) {
  lv <- levels(g); a <- x[g == lv[1]]; b <- x[g == lv[2]]
  est <- function(a, b) {
    r <- rank(c(a, b))
    u <- sum(r[length(a) + seq_along(b)]) - length(b) * (length(b) + 1) / 2
    2 * u / (length(a) * length(b)) - 1
  }
  bs <- replicate(1000, est(sample(a, replace = TRUE), sample(b, replace = TRUE)))
  c(est(a, b), unname(quantile(bs, c(0.025, 0.975))))
}

eps2 <- function(x, g) unname(kruskal.test(x ~ g)$statistic) / (length(x) - 1)
eps2_ci <- function(x, g) {
  b <- replicate(n_boot, {
    i <- sample(length(x), replace = TRUE); gi <- droplevels(g[i])
    if (nlevels(gi) < 2) NA_real_ else eps2(x[i], gi)
  })
  unname(quantile(b, c(0.025, 0.975), na.rm = TRUE))
}

spearman_ci <- function(x, y) {
  n <- length(x); r <- cor(x, y, method = "spearman")
  se <- sqrt((1 + r^2 / 2) / (n - 3)); z <- atanh(r)
  c(r, tanh(z - 1.96 * se), tanh(z + 1.96 * se))
}

mk <- function(base, comparison = NA, n = NA, gn = NA, test = NA, em = NA,
               e = c(NA, NA, NA), p = NA, note = "") {
  as.data.frame(c(base, list(comparison = comparison, n = n, group_n = gn,
    test = test, effect_measure = em, effect = e[1], ci_low = e[2], ci_high = e[3],
    p_value = p, note = note)), stringsAsFactors = FALSE)
}

test_one <- function(cancer, dat, f, omics, v, type) {
  base <- list(cancer = cancer, omics = omics, feature = f, variable = v)
  x <- dat[[f]]; y <- dat[[v]]
  if (type == "cont") {
    ok <- !is.na(x) & !is.na(y)
    if (sum(ok) < min_cont_n) return(mk(base, n = sum(ok), note = paste0("not tested: n<", min_cont_n)))
    p <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)$p.value
    return(mk(base, "correlation with age (years)", sum(ok), NA, "Spearman",
              "rho", spearman_ci(x[ok], y[ok]), p))
  }
  ok <- !is.na(x) & !is.na(y)
  g <- droplevels(factor(y[ok])); xx <- x[ok]
  tab <- table(g)
  gn <- paste(names(tab), tab, sep = "=", collapse = "; ")
  if (type == "two") {
    if (nlevels(g) < 2 || min(tab) < min_group_n)
      return(mk(base, n = sum(ok), gn = gn, note = paste0("not tested: smallest group n<", min_group_n)))
    p <- wilcox.test(xx ~ g, exact = FALSE)$p.value
    return(mk(base, paste(levels(g)[2], "vs", levels(g)[1]), sum(ok), gn,
              "Wilcoxon rank-sum", "Cliff delta", cliff_delta(xx, g), p))
  }
  keep <- names(tab)[tab >= min_group_n]
  note <- if (length(keep) < length(tab)) "groups with n<5 excluded" else ""
  sel <- g %in% keep; g <- droplevels(g[sel]); xx <- xx[sel]
  if (nlevels(g) < 2) return(mk(base, n = sum(ok), gn = gn, note = "not tested: <2 usable groups"))
  p <- kruskal.test(xx ~ g)$p.value
  ci <- eps2_ci(xx, g)
  mk(base, paste(nlevels(g), "groups"), length(xx), gn, "Kruskal-Wallis",
     "epsilon-squared", c(eps2(xx, g), ci[1], ci[2]), p, note)
}

var_specs <- list(
  BRCA = list(age = "cont", stage = "two", nodal = "two", subtype = "multi"),
  PRAD = list(age = "cont", gleason = "multi", highgrade = "two", nodal = "two"))

data_store <- list(); results <- list()
for (cn in c("BRCA", "PRAD")) {
  ds <- build_dataset(cn); data_store[[cn]] <- ds
  for (i in seq_len(nrow(ds$features))) {
    for (v in names(var_specs[[cn]])) {
      results[[length(results) + 1]] <- test_one(cn, ds$data, ds$features$feature[i],
        ds$features$omics[i], v, var_specs[[cn]][[v]])
    }
  }
}
res <- do.call(rbind, results)
res$fdr <- NA_real_
for (cn in unique(res$cancer)) {
  i <- res$cancer == cn & !is.na(res$p_value)
  res$fdr[i] <- p.adjust(res$p_value[i], "BH")
}
write.csv(res, file.path(out_dir, "Step9_1_clinical_association_results.csv"), row.names = FALSE)

# missing-data documentation
miss <- do.call(rbind, lapply(names(data_store), function(cn) {
  d <- data_store[[cn]]$clinical; vars <- setdiff(names(d), "patient")
  data.frame(cancer = cn, variable = vars, n_patients = nrow(d),
             n_available = vapply(vars, function(v) sum(!is.na(d[[v]])), numeric(1)),
             stringsAsFactors = FALSE)
}))
miss$pct_missing <- round(100 * (1 - miss$n_available / miss$n_patients), 1)
extra <- data.frame(cancer = c("BRCA", "PRAD"), variable = c("tumor grade", "pathologic stage"),
                    n_patients = c(nrow(data_store$BRCA$clinical), nrow(data_store$PRAD$clinical)),
                    n_available = 0, pct_missing = 100)
miss <- rbind(miss, extra)
write.csv(miss, file.path(out_dir, "Step9_1_missing_clinical_data.csv"), row.names = FALSE)

feat_n <- do.call(rbind, lapply(names(data_store), function(cn) {
  d <- data_store[[cn]]$data; f <- data_store[[cn]]$features
  data.frame(cancer = cn, feature = f$feature, omics = f$omics,
             n_with_value = vapply(f$feature, function(z) sum(!is.na(d[[z]])), numeric(1)))
}))
write.csv(feat_n, file.path(out_dir, "Step9_1_feature_availability.csv"), row.names = FALSE)

# plots for FDR < 0.05
for (cn in names(data_store)) {
  s <- res[res$cancer == cn & !is.na(res$fdr) & res$fdr < 0.05, ]
  if (!nrow(s)) next
  d <- data_store[[cn]]$data
  pdf(file.path(out_dir, paste0("Step9_1_", cn, "_significant.pdf")), width = 5, height = 5)
  for (i in seq_len(nrow(s))) {
    f <- s$feature[i]; v <- s$variable[i]
    ttl <- sprintf("%s ~ %s (FDR=%.3g)", f, v, s$fdr[i])
    if (is.numeric(d[[v]])) plot(d[[v]], d[[f]], xlab = v, ylab = f, main = ttl, pch = 19, cex = 0.5)
    else boxplot(d[[f]] ~ d[[v]], xlab = "", ylab = f, main = ttl, las = 2)
  }
  dev.off()
}

cat("\n=== STEP 9.1 SUMMARY ===\n")
cat("Tests run:", sum(!is.na(res$p_value)), "| not tested:", sum(is.na(res$p_value)),
    "| FDR<0.05:", sum(res$fdr < 0.05, na.rm = TRUE), "\n")
cat("\nNot tested:\n")
print(res[is.na(res$p_value), c("cancer", "feature", "variable", "note")], row.names = FALSE)
cat("\nTwo-group tests (stage / nodal / Gleason 6-7 vs 8-9), sorted by p:\n")
two <- res[res$effect_measure %in% "Cliff delta", c("cancer", "feature", "variable", "n",
            "effect", "ci_low", "ci_high", "p_value", "fdr")]
two <- two[order(two$p_value), ]
two[sapply(two, is.numeric)] <- lapply(two[sapply(two, is.numeric)], function(z) signif(z, 3))
print(head(two, 15), row.names = FALSE)
cat("\nSaved in:", out_dir, "\n")
