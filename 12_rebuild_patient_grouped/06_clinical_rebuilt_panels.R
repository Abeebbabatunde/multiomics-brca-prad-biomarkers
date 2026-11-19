# 06_clinical_rebuilt_panels.R
# Step 12.6 (Objective 2, step I): clinical association (9.1) and overall-survival (9.2) analyses of the thesis, re-run UNCHANGED on the REBUILT panels.
# The two thesis scripts (09_clinical/09_1_clinical_association_RPPA_EDITED.R and 09_2_survival.R) are included below with three edits only:
#   (1) the panel is the rebuilt panel (BRCA: 2 CpGs + 10 genes; PRAD: 3 genes; PRAD methylation panel is empty), genes mapped by their exact Ensembl IDs;
#   (2) the protein feature (VAV1) is not part of the rebuild;  (3) outputs go to 12_rebuild_clinical/.
# Same data, tests, thresholds, FDR control and events-per-covariate rules as the thesis. Gene expression for the clinical analysis comes from the thesis
# processed files (all tumours), as in the thesis; the rebuild transformation (vst blind = TRUE) applies to the 56 / 17 pairs only.
# With 10 breast genes, the multivariable "RNA panel vs clinical" model needs >= 130 events and will be reported as not fitted if events are fewer.
# Edit project_root only (inside the first lines below).

# STEP 9.1 CLINICAL ASSOCIATION - PART 1/3
# EDITED COPY: only get_protein() differs from the run version (reads RPPA instead of the iTRAQ file). Everything else is unchanged.
project_root <- "C:/Users/USER/Documents/MyThesis_MultiOmics"
P <- function(...) file.path(project_root, ...)
out_dir <- P("12_rebuild_clinical")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(42)
min_group_n <- 5      # smallest group allowed in a comparison
min_cont_n  <- 20     # smallest n for a correlation
n_boot      <- 1000   # bootstrap resamples for epsilon-squared CIs
## ---- REBUILT panels (replaces the locked Step 7 panel) ----
rbp <- function(nm) read.csv(file.path(project_root, paste0(nm, "_rebuilt_panel.csv")), stringsAsFactors = FALSE)
b_m <- rbp("BRCA_methylation"); b_r <- rbp("BRCA_RNA"); p_r <- rbp("PRAD_RNA")
panel_rb <- list(BRCA = c(b_m$feature, b_r$symbol), PRAD = p_r$symbol)          # PRAD methylation: empty rebuilt panel
# genes are mapped to Ensembl by the exact IDs in the rebuilt panel (no symbol / alias look-up)
manual_ensembl <- c(setNames(sub("\\..*$", "", b_r$feature), b_r$symbol), setNames(sub("\\..*$", "", p_r$feature), p_r$symbol))

for (pk in c("data.table", "org.Hs.eg.db", "AnnotationDbi",
             "SummarizedExperiment", "DESeq2")) {
  if (!requireNamespace(pk, quietly = TRUE))
    stop("Install this package first: ", pk, call. = FALSE)
}

panel <- panel_rb
protein_feats <- list(BRCA = character(0), PRAD = character(0))   # protein layer is not part of the rebuild

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


## =================== survival (thesis 9.2) ===================
# STEP 9.2 SURVIVAL (BRCA overall survival) - PART 1/2
suppressPackageStartupMessages(library(survival))
out_dir <- P("12_rebuild_clinical")
set.seed(42)
EPV <- 10   # events required per model covariate

os_table <- function(cancer) {
  f <- if (cancer == "BRCA") P("01_RawData", "Breast", "clinical_brca.csv")
       else P("01_RawData", "Prostate", "clinical_prad.csv")
  cl <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  dead <- as.integer(tolower(cl$vital_status) == "dead")
  t <- suppressWarnings(as.numeric(ifelse(dead == 1, cl$days_to_death, cl$days_to_last_follow_up)))
  data.frame(patient = cl$submitter_id, os_days = t, os_event = dead, stringsAsFactors = FALSE)
}
usable_os <- function(os) os[!is.na(os$os_days) & !is.na(os$os_event) & os$os_days > 0, ]

# ---- endpoint audit (both cancers)
audit <- do.call(rbind, lapply(c("BRCA", "PRAD"), function(cn) {
  os <- os_table(cn); d <- usable_os(os)
  rv <- survfit(Surv(os_days, 1 - os_event) ~ 1, data = d)   # reverse KM follow-up
  data.frame(cancer = cn, n_clinical = nrow(os), n_usable = nrow(d),
    events = sum(d$os_event), pct_events = round(100 * mean(d$os_event), 1),
    median_followup_yrs = round(unname(summary(rv)$table["median"]) / 365.25, 2),
    max_followup_yrs = round(max(d$os_days) / 365.25, 1),
    max_covariates_10EPV = floor(sum(d$os_event) / 10))
}))
write.csv(audit, file.path(out_dir, "Step9_2_endpoint_audit.csv"), row.names = FALSE)
cat("\n=== OVERALL SURVIVAL AUDIT ===\n"); print(audit, row.names = FALSE)
cat("\nPRAD: too few events for any Cox model; reported as not feasible.\n")

# ---- BRCA analysis dataset
ds <- if (exists("data_store") && !is.null(data_store$BRCA)) data_store$BRCA else build_dataset("BRCA")
dat <- merge(ds$data, usable_os(os_table("BRCA")), by = "patient")
dat$os_years <- dat$os_days / 365.25
dat$age10 <- dat$age / 10
dat$stage <- relevel(factor(dat$stage), ref = "Stage I")
dat$nodal <- relevel(factor(dat$nodal), ref = "N0")
feats <- ds$features[ds$features$omics != "protein", ]   # protein features excluded (none in the rebuilt panels)
zs <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
for (f in feats$feature) dat[[paste0("z_", f)]] <- zs(dat[[f]])
cat("\nBRCA survival cohort:", nrow(dat), "patients,", sum(dat$os_event), "deaths\n")

sf <- function(zc) paste("Surv(os_years, os_event) ~", zc)
ph_last <- function(fit) tryCatch({ t <- cox.zph(fit)$table; unname(t[nrow(t), "p"]) },
                                  error = function(e) NA_real_)

# ---- univariate Cox (HR per 1 SD)
uni <- do.call(rbind, lapply(seq_len(nrow(feats)), function(i) {
  f <- feats$feature[i]; zc <- paste0("z_", f)
  d <- dat[!is.na(dat[[zc]]), ]
  fit <- coxph(as.formula(sf(zc)), data = d); s <- summary(fit)
  data.frame(feature = f, omics = feats$omics[i], n = nrow(d), events = sum(d$os_event),
    HR_per_SD = s$conf.int[1, 1], ci_low = s$conf.int[1, 3], ci_high = s$conf.int[1, 4],
    p_value = s$coefficients[1, 5], ph_global_p = ph_last(fit))
}))
uni$fdr <- p.adjust(uni$p_value, "BH")
write.csv(uni, file.path(out_dir, "Step9_2_cox_univariate.csv"), row.names = FALSE)

# ---- multivariable Cox: biomarker + age (per 10 y) + stage + nodal
adj <- do.call(rbind, lapply(seq_len(nrow(feats)), function(i) {
  f <- feats$feature[i]; zc <- paste0("z_", f)
  d <- dat[complete.cases(dat[, c(zc, "age10", "stage", "nodal")]), ]
  ev <- sum(d$os_event)
  out <- data.frame(feature = f, omics = feats$omics[i], n = nrow(d), events = ev,
    HR_per_SD = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_,
    ph_global_p = NA_real_, ph_biomarker_p = NA_real_, note = "")
  if (ev < EPV * 4) { out$note <- "not fitted: <10 events per covariate"; return(out) }
  fit <- coxph(as.formula(paste(sf(zc), "+ age10 + stage + nodal")), data = d)
  s <- summary(fit); zp <- tryCatch(cox.zph(fit)$table, error = function(e) NULL)
  out$HR_per_SD <- s$conf.int[zc, 1]; out$ci_low <- s$conf.int[zc, 3]
  out$ci_high <- s$conf.int[zc, 4]; out$p_value <- s$coefficients[zc, 5]
  if (!is.null(zp)) { out$ph_global_p <- unname(zp["GLOBAL", "p"]); out$ph_biomarker_p <- unname(zp[zc, "p"]) }
  out
}))
adj$fdr <- p.adjust(adj$p_value, "BH")
write.csv(adj, file.path(out_dir, "Step9_2_cox_multivariable.csv"), row.names = FALSE)
cat("Cox models done. Now run Part 2.\n")

# PART 2/2

# ---- Kaplan-Meier, median split (fixed in advance), all features plotted
km_rows <- list()
pdf(file.path(out_dir, "Step9_2_KM_median_split.pdf"), width = 6, height = 5)
for (i in seq_len(nrow(feats))) {
  f <- feats$feature[i]
  d <- dat[!is.na(dat[[f]]), ]
  d$grp <- factor(ifelse(d[[f]] > median(d[[f]]), "High", "Low"), levels = c("Low", "High"))
  ev <- tapply(d$os_event, d$grp, sum)
  row <- data.frame(feature = f, omics = feats$omics[i], n = nrow(d),
    events_low = as.integer(ev["Low"]), events_high = as.integer(ev["High"]),
    p_logrank = NA_real_, note = "", stringsAsFactors = FALSE)
  if (any(is.na(ev)) || min(ev) < 10) {
    row$note <- "not run: <10 events in a group"
  } else {
    p <- 1 - pchisq(survdiff(Surv(os_years, os_event) ~ grp, data = d)$chisq, 1)
    row$p_logrank <- p
    fit <- survfit(Surv(os_years, os_event) ~ grp, data = d)
    plot(fit, col = c("steelblue", "firebrick"), lwd = 2, mark.time = TRUE,
         xlab = "Years", ylab = "Overall survival",
         main = sprintf("%s (median split)\nlog-rank p = %.3g", f, p))
    legend("bottomleft", bty = "n", lwd = 2, col = c("steelblue", "firebrick"),
      legend = sprintf("%s (n=%d, events=%d)", levels(d$grp), as.integer(table(d$grp)), as.integer(ev)))
  }
  km_rows[[i]] <- row
}
dev.off()
km <- do.call(rbind, km_rows)
km$fdr <- p.adjust(km$p_logrank, "BH")
write.csv(km, file.path(out_dir, "Step9_2_KM_logrank.csv"), row.names = FALSE)

# ---- RNA panel model: does the 3-gene RNA panel add to age + stage + nodal?
rna <- feats$feature[feats$omics == "RNA"]
zr <- paste0("z_", rna)
d <- dat[complete.cases(dat[, c(zr, "age10", "stage", "nodal")]), ]
ev <- sum(d$os_event); ncov <- length(zr) + 3
panel_out <- NULL
if (ev >= EPV * ncov) {
  clin <- coxph(Surv(os_years, os_event) ~ age10 + stage + nodal, data = d)
  full <- coxph(as.formula(paste("Surv(os_years, os_event) ~ age10 + stage + nodal +",
                                 paste(zr, collapse = " + "))), data = d)
  lr <- anova(clin, full)
  sm <- summary(full)
  terms_tbl <- data.frame(term = rownames(sm$conf.int), HR = sm$conf.int[, 1],
    ci_low = sm$conf.int[, 3], ci_high = sm$conf.int[, 4],
    p_value = sm$coefficients[, 5], row.names = NULL)
  write.csv(terms_tbl, file.path(out_dir, "Step9_2_RNA_panel_cox_terms.csv"), row.names = FALSE)
  zp <- tryCatch(cox.zph(full)$table, error = function(e) NULL)
  panel_out <- data.frame(n = nrow(d), events = ev,
    C_clinical_only = summary(clin)$concordance[1], C_clinical_plus_RNA = sm$concordance[1],
    LR_p_value = lr[2, "Pr(>|Chi|)"],
    ph_global_p = if (is.null(zp)) NA else unname(zp["GLOBAL", "p"]))
  write.csv(panel_out, file.path(out_dir, "Step9_2_RNA_panel_vs_clinical.csv"), row.names = FALSE)
} else cat("RNA panel model not fitted: events < 10 per covariate\n")

# ---- console summary
rd <- function(df) { n <- sapply(df, is.numeric); df[n] <- lapply(df[n], signif, 3); df }
cat("\n=== UNIVARIATE COX (HR per 1 SD) ===\n"); print(rd(uni), row.names = FALSE)
cat("\n=== ADJUSTED COX (age, stage, nodal) ===\n"); print(rd(adj), row.names = FALSE)
cat("\n=== KAPLAN-MEIER LOG-RANK ===\n"); print(rd(km), row.names = FALSE)
if (!is.null(panel_out)) { cat("\n=== RNA PANEL vs CLINICAL-ONLY ===\n"); print(rd(panel_out), row.names = FALSE) }
cat("\nSaved in:", out_dir, "\n")
