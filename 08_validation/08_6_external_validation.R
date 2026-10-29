# STEP 8.6-8.9 EXTERNAL VALIDATION (robust version) - PART 1/3
project_root <- "C:/Users/USER/Documents/MyThesis_MultiOmics"

cfg8e <- list(
  step7_dir     = file.path(project_root, "07_FinalResults", "Step7_FeatureSelection_ML"),
  step8_dir     = file.path(project_root, "08_FinalResults", "Step8_Validation"),
  geo_cache_dir = file.path(project_root, "08_FinalResults", "External_GEO_cache"),
  out_dir       = file.path(project_root, "08_FinalResults", "Step8_External_Validation"),
  seed     = 42,
  rf_ntree = 1000
)
dir.create(cfg8e$out_dir,       recursive = TRUE, showWarnings = FALSE)
dir.create(cfg8e$geo_cache_dir, recursive = TRUE, showWarnings = FALSE)

options(timeout = max(600, getOption("timeout")))

suppressPackageStartupMessages({
  library(Biobase)
  library(GEOquery)
  library(randomForest)
  library(caret)
  library(pROC)
})

require_file <- function(path) {
  if (!file.exists(path)) stop("Required file not found: ", path, call. = FALSE)
  path
}

with_retry <- function(fun, what, tries = 4, wait = 5) {
  last_err <- NULL
  for (i in seq_len(tries)) {
    res <- tryCatch(list(ok = TRUE, value = fun()),
                    error = function(e) list(ok = FALSE, err = conditionMessage(e)))
    if (res$ok) return(res$value)
    last_err <- res$err
    message(sprintf("[retry %d/%d] %s failed: %s", i, tries, what, last_err))
    if (i < tries) Sys.sleep(wait * i)
  }
  stop(sprintf("%s failed after %d attempts. Last error: %s", what, tries, last_err),
       call. = FALSE)
}

geo_stub <- function(id) sub("\\d{1,3}$", "nnn", id)

download_robust <- function(url, dest) {
  tmp <- paste0(dest, ".part")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  good <- function() file.exists(tmp) && isTRUE(file.info(tmp)$size > 0)
  ok <- FALSE
  if (requireNamespace("curl", quietly = TRUE)) {
    ok <- tryCatch({
      h <- curl::new_handle(http_version = 2L, timeout = 600, connecttimeout = 60)
      curl::curl_download(url, tmp, handle = h, quiet = TRUE)
      good()
    }, error = function(e) {
      message("  curl (HTTP/1.1) failed: ", conditionMessage(e)); FALSE
    })
  }
  if (!ok) {
    ok <- tryCatch({
      suppressWarnings(utils::download.file(url, tmp, mode = "wb", quiet = TRUE, method = "libcurl"))
      good()
    }, error = function(e) {
      message("  download.file failed: ", conditionMessage(e)); FALSE
    })
  }
  if (!ok) return(FALSE)
  isTRUE(file.copy(tmp, dest, overwrite = TRUE))
}

load_series_matrix <- function(gse_id) {
  f <- file.path(cfg8e$geo_cache_dir, paste0(gse_id, "_series_matrix.txt.gz"))
  fetch <- function() {
    url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s_series_matrix.txt.gz",
                   geo_stub(gse_id), gse_id, gse_id)
    cat("Downloading", gse_id, "series matrix...\n")
    with_retry(function() if (!download_robust(url, f)) stop("download failed"),
               what = paste("Download of", gse_id), tries = 3)
  }
  if (!file.exists(f) || isTRUE(file.info(f)$size == 0)) fetch()
  es <- tryCatch(
    GEOquery::getGEO(filename = f, getGPL = FALSE),
    error = function(e) {
      message("Cached file for ", gse_id, " unreadable (", conditionMessage(e),
              "); re-downloading once.")
      unlink(f); fetch()
      GEOquery::getGEO(filename = f, getGPL = FALSE)
    }
  )
  if (!inherits(es, "ExpressionSet") && is.list(es)) es <- es[[1]]
  es
}

parse_platform_file <- function(f) {
  gpl <- GEOquery::getGEO(filename = f)
  tbl <- GEOquery::Table(gpl)
  if (is.null(tbl) || nrow(tbl) == 0) stop("platform table is empty")
  as.data.frame(tbl, stringsAsFactors = FALSE)
}

get_platform_table <- function(gpl_id) {
  cache <- cfg8e$geo_cache_dir
  local_files <- file.path(cache, c(paste0(gpl_id, "_family.soft.gz"),
                                    paste0(gpl_id, ".soft.gz"),
                                    paste0(gpl_id, ".soft")))
  local_files <- local_files[file.exists(local_files)]
  for (f in local_files) {
    tbl <- tryCatch(parse_platform_file(f), error = function(e) {
      message("Cached platform file unreadable (", basename(f), "): ", conditionMessage(e))
      NULL
    })
    if (!is.null(tbl)) { cat("Using cached platform file:", basename(f), "\n"); return(tbl) }
  }
  dest <- file.path(cache, paste0(gpl_id, "_family.soft.gz"))
  url  <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/platforms/%s/%s/soft/%s_family.soft.gz",
                  geo_stub(gpl_id), gpl_id, gpl_id)
  cat("Downloading platform annotation for", gpl_id, "(HTTP/1.1)...\n")
  got <- tryCatch(
    with_retry(function() if (!download_robust(url, dest)) stop("download failed"),
               what = paste("Download of", gpl_id), tries = 3),
    error = function(e) { message(conditionMessage(e)); NULL })
  if (!is.null(got)) {
    tbl <- tryCatch(parse_platform_file(dest), error = function(e) {
      message("Downloaded platform file unreadable: ", conditionMessage(e)); unlink(dest); NULL
    })
    if (!is.null(tbl)) return(tbl)
  }
  tryCatch(
    with_retry(function() as.data.frame(GEOquery::Table(GEOquery::getGEO(gpl_id, destdir = cache))),
               what = paste("GEOquery fetch of", gpl_id), tries = 3),
    error = function(e) stop(
      "Could not obtain platform annotation ", gpl_id, ".\n",
      "MANUAL FIX: download\n  ", url, "\nand save it as\n  ", dest, "\nthen re-run.\n",
      "(Original error: ", conditionMessage(e), ")", call. = FALSE)
  )
}

# PART 2/3

map_probes_to_symbols <- function(eset) {
  gpl_id <- Biobase::annotation(eset)
  if (is.null(gpl_id) || !nzchar(gpl_id)) stop("No platform ID found in the series matrix header.")
  tbl <- get_platform_table(gpl_id)
  nm  <- names(tbl)

  pick <- function(patterns) {
    for (p in patterns) {
      hit <- nm[grepl(p, nm, ignore.case = TRUE)]
      if (length(hit)) return(hit[1])
    }
    NA_character_
  }

  symbols <- NULL
  sym_col <- pick(c("^gene[ _.]?symbol$", "^symbol$", "^ilmn_gene$", "gene[ _.]?symbol"))
  if (!is.na(sym_col)) {
    symbols <- as.character(tbl[[sym_col]])
  } else {
    assign_col <- pick(c("gene_assignment", "gene assignment"))
    if (!is.na(assign_col)) {
      extract_symbols <- function(s) {
        if (is.na(s) || s == "" || s == "---") return(NA_character_)
        entries <- strsplit(s, "///", fixed = TRUE)[[1]]
        syms <- vapply(entries, function(e) {
          toks <- trimws(strsplit(e, "//", fixed = TRUE)[[1]])
          if (length(toks) >= 2) toks[2] else NA_character_
        }, character(1))
        syms <- syms[!is.na(syms)]
        if (length(syms) == 0) NA_character_ else paste(unique(syms), collapse = " /// ")
      }
      symbols <- vapply(tbl[[assign_col]], extract_symbols, character(1), USE.NAMES = FALSE)
    }
  }
  if (is.null(symbols)) {
    stop("No gene symbol column found in platform ", gpl_id,
         ". Columns are: ", paste(nm, collapse = ", "), call. = FALSE)
  }

  id_col <- if ("ID" %in% nm) "ID" else nm[1]
  raw <- data.frame(probe_id = as.character(tbl[[id_col]]), symbol = symbols,
                    stringsAsFactors = FALSE)
  raw <- raw[!is.na(raw$symbol) & nzchar(trimws(raw$symbol)) & raw$symbol != "---", , drop = FALSE]

  parts <- strsplit(raw$symbol, "\\s*///\\s*|\\s*;\\s*")
  long <- data.frame(probe_id = rep(raw$probe_id, lengths(parts)),
                     symbol   = trimws(unlist(parts)),
                     stringsAsFactors = FALSE)
  long <- unique(long[nzchar(long$symbol) & long$symbol != "---", , drop = FALSE])
  cat("Mapped", nrow(long), "probe-gene pairs from", gpl_id, "\n")
  long
}

collapse_to_panel_genes <- function(eset, genes) {
  expr <- Biobase::exprs(eset)
  storage.mode(expr) <- "double"

  qx <- stats::quantile(expr, c(0, 0.99), na.rm = TRUE)
  if (qx[2] > 100) {
    message("Expression values look un-logged (99th pct = ", round(qx[2], 1),
            "); applying log2(x + 1).")
    expr <- log2(pmax(expr, 0) + 1)
  }

  map <- map_probes_to_symbols(eset)
  map <- map[map$symbol %in% genes & map$probe_id %in% rownames(expr), , drop = FALSE]
  if (nrow(map) == 0) stop("None of the panel genes could be mapped to probes on this platform.")

  probe_var <- apply(expr[unique(map$probe_id), , drop = FALSE], 1, stats::var, na.rm = TRUE)
  map$v <- probe_var[map$probe_id]
  map <- map[order(-map$v), , drop = FALSE]
  best <- map[!duplicated(map$symbol), , drop = FALSE]

  out <- t(expr[best$probe_id, , drop = FALSE])
  colnames(out) <- best$symbol
  as.data.frame(out)
}

extract_meth_matrix <- function(eset, feats) {
  m <- Biobase::exprs(eset)
  present <- intersect(feats, rownames(m))
  if (length(present) == 0) stop("None of the panel CpG probes are present in this dataset.")
  sub <- m[present, , drop = FALSE]
  rng <- range(sub, na.rm = TRUE)
  if (rng[1] < 0 || rng[2] > 1) {
    message("WARNING: methylation values outside [0,1] (range ", round(rng[1], 2), " to ",
            round(rng[2], 2), ") - these may be M-values, not betas.")
  }
  as.data.frame(t(sub))
}

get_labels <- function(pdata, col, mapping, gse_id) {
  if (!col %in% names(pdata)) {
    cand <- grep(":ch1$", names(pdata), value = TRUE)
    stop(sprintf("Column '%s' not found in %s. Available characteristic columns: %s",
                 col, gse_id, paste(cand, collapse = ", ")), call. = FALSE)
  }
  raw <- trimws(as.character(pdata[[col]]))
  names(mapping) <- tolower(names(mapping))
  lab <- unname(mapping[tolower(raw)])
  bad <- unique(raw[is.na(lab)])
  if (length(bad) > 0) {
    message(gse_id, ": dropping samples with unmapped labels: ", paste(bad, collapse = " | "))
  }
  stats::setNames(lab, rownames(pdata))
}

zscore_cols <- function(X) {
  X <- as.matrix(X)
  m <- colMeans(X, na.rm = TRUE)
  s <- apply(X, 2, stats::sd, na.rm = TRUE)
  s[is.na(s) | s == 0] <- 1
  as.data.frame(sweep(sweep(X, 2, m, "-"), 2, s, "/"))
}

train_single_omics_model <- function(full_df, feats, train_patients, zscore = FALSE, seed) {
  miss <- setdiff(feats, names(full_df))
  if (length(miss) > 0) stop("Features missing from internal dataset: ", paste(miss, collapse = ", "))
  train_mask <- full_df$patient %in% train_patients
  if (!any(train_mask)) stop("No internal training patients matched the locked train list.")

  X <- as.data.frame(full_df[train_mask, feats, drop = FALSE])
  X[] <- lapply(X, as.numeric)
  y <- factor(full_df$group[train_mask], levels = c("Normal", "Early"))

  ok <- stats::complete.cases(X) & !is.na(y)
  if (any(!ok)) {
    message("Training: dropping ", sum(!ok), " internal sample(s) with NA values/labels.")
    X <- X[ok, , drop = FALSE]; y <- y[ok]
  }
  if (nlevels(droplevels(y)) < 2) stop("Internal training data has only one class.")

  if (zscore) X <- zscore_cols(X)

  set.seed(seed)
  model <- randomForest::randomForest(x = X, y = y, ntree = cfg8e$rf_ntree)
  list(model = model, feats = feats, zscore = zscore)
}

apply_single_omics_model <- function(trained, X_external) {
  X <- as.data.frame(X_external[, trained$feats, drop = FALSE])
  X[] <- lapply(X, as.numeric)
  if (trained$zscore) X <- zscore_cols(X)
  prob <- predict(trained$model, X, type = "prob")[, "Early"]
  pred <- predict(trained$model, X, type = "response")
  list(pred = pred, prob = prob)
}

evaluate_external <- function(pred, prob, true_labels, sample_ids, cohort_name, omics_name) {
  true_f <- factor(true_labels, levels = c("Normal", "Early"))
  pred_f <- factor(as.character(pred), levels = c("Normal", "Early"))
  cm <- caret::confusionMatrix(pred_f, true_f, positive = "Early")
  auc_val <- tryCatch(
    as.numeric(pROC::auc(pROC::roc(response = true_f, predictor = as.numeric(prob),
                                   levels = c("Normal", "Early"), direction = "<", quiet = TRUE))),
    error = function(e) NA_real_)

  cat("\n---", cohort_name, "-", omics_name, "external performance (n =", length(true_f), ") ---\n")
  print(cm)
  cat("AUC:", auc_val, "\n")

  write.csv(data.frame(sample = sample_ids, true = as.character(true_f),
                       pred = as.character(pred_f), prob_early = as.numeric(prob)),
            file.path(cfg8e$out_dir, paste0(cohort_name, "_", omics_name, "_external_predictions.csv")),
            row.names = FALSE)

  list(confusion_matrix = cm, auc = auc_val, n = length(true_f))
}

run_external <- function(cohort, omics, locked_feats, full_df, train_patients,
                         ext_X, ext_y, zscore) {
  ext_y <- ext_y[rownames(ext_X)]
  keep  <- !is.na(ext_y)
  ext_X <- ext_X[keep, , drop = FALSE]
  ext_y <- ext_y[keep]
  cat("External label counts:", paste(names(table(ext_y)), table(ext_y), collapse = " / "), "\n")

  feats_used <- intersect(locked_feats, colnames(ext_X))
  missing    <- setdiff(locked_feats, feats_used)
  if (length(missing) > 0) {
    cat("WARNING: locked features missing in external data:", paste(missing, collapse = ", "),
        "\n  -> model trained on the", length(feats_used), "shared features only (disclose in Methods).\n")
  }
  if (length(feats_used) == 0) stop("No locked features available in the external data.")

  ext_X <- as.data.frame(ext_X[, feats_used, drop = FALSE])
  ext_X[] <- lapply(ext_X, as.numeric)
  ok <- stats::complete.cases(ext_X)
  if (any(!ok)) {
    message("Dropping ", sum(!ok), " external sample(s) with NA in panel features.")
    ext_X <- ext_X[ok, , drop = FALSE]; ext_y <- ext_y[ok]
  }
  if (length(unique(ext_y)) < 2) stop("External cohort has only one class after filtering.")

  trained <- train_single_omics_model(full_df, feats_used, train_patients, zscore, cfg8e$seed)
  pr  <- apply_single_omics_model(trained, ext_X)
  res <- evaluate_external(pr$pred, pr$prob, ext_y, rownames(ext_X), cohort, omics)
  res$features_used    <- feats_used
  res$features_missing <- missing
  res
}

safe_run <- function(label, fun) {
  tryCatch(fun(), error = function(e) {
    cat("\n!! ", label, " FAILED: ", conditionMessage(e), "\n", sep = "")
    list(error = conditionMessage(e), auc = NA_real_, n = NA_integer_)
  })
}

banner <- function(txt) {
  cat("\n============================================================\n")
  cat(txt, "\n")
  cat("============================================================\n")
}

# PART 3/3

step7_results <- readRDS(require_file(file.path(cfg8e$step7_dir, "Step7_feature_selection_results.rds")))
brca_full <- readRDS(require_file(file.path(cfg8e$step7_dir, "BRCA_targeted_dataset.rds")))
prad_full <- readRDS(require_file(file.path(cfg8e$step7_dir, "PRAD_targeted_dataset.rds")))

brca_train_patients <- read.csv(require_file(file.path(cfg8e$step8_dir, "BRCA_locked_train_patients.csv")))$patient
prad_train_patients <- read.csv(require_file(file.path(cfg8e$step8_dir, "PRAD_locked_train_patients.csv")))$patient

brca_panel <- step7_results$BRCA$consensus
prad_panel <- step7_results$PRAD$consensus
if (length(brca_panel) == 0 || length(prad_panel) == 0) stop("Locked panel is empty - check Step 7 results.")

brca_consensus_genes <- read.csv(require_file(file.path(project_root, "06_FinalResults", "Step6_Network",
                                                         "BRCA_3omics", "consensus_genes.csv")),
                                 stringsAsFactors = FALSE)

block_of <- function(feature, consensus_df) {
  if (!"block" %in% names(consensus_df)) stop("'block' column missing from consensus_genes.csv")
  idx <- rep(FALSE, nrow(consensus_df))
  for (col in intersect(c("gene", "feature"), names(consensus_df))) {
    idx <- idx | ((as.character(consensus_df[[col]]) == feature) %in% TRUE)
  }
  hit <- consensus_df$block[idx]
  if (length(hit) == 0) NA_character_ else tolower(as.character(hit[1]))
}

brca_meth_feats  <- brca_panel[grepl("^cg", brca_panel)]
brca_other_feats <- setdiff(brca_panel, brca_meth_feats)
brca_blocks      <- vapply(brca_other_feats, block_of, character(1), consensus_df = brca_consensus_genes)
if (anyNA(brca_blocks)) {
  warning("No block found for BRCA feature(s): ",
          paste(brca_other_feats[is.na(brca_blocks)], collapse = ", "))
}
brca_rna_feats  <- brca_other_feats[brca_blocks %in% "rna"]
brca_prot_feats <- brca_other_feats[brca_blocks %in% c("prot", "protein")]

prad_meth_feats <- prad_panel[grepl("^cg", prad_panel)]
prad_rna_feats  <- setdiff(prad_panel, prad_meth_feats)

cat("BRCA panel split -- meth:", length(brca_meth_feats), "| rna:", length(brca_rna_feats),
    "| prot (excluded from external):", length(brca_prot_feats), "\n")
cat("PRAD panel split -- meth:", length(prad_meth_feats), "| rna:", length(prad_rna_feats), "\n")

# ---- BRCA methylation: GSE60185
brca_meth_ext_result <- safe_run("BRCA methylation (GSE60185)", function() {
  banner("BRCA METHYLATION EXTERNAL VALIDATION - GSE60185")
  es <- load_series_matrix("GSE60185")
  pd <- Biobase::pData(es)
  y  <- get_labels(pd, "disease state:ch1",
                   c("normal" = "Normal", "breast cancer" = "Early"), "GSE60185")
  X  <- extract_meth_matrix(es, brca_meth_feats)
  run_external("BRCA", "methylation", brca_meth_feats, brca_full, brca_train_patients,
               X, y, zscore = FALSE)
})

# ---- BRCA RNA: GSE54002
brca_rna_ext_result <- safe_run("BRCA RNA (GSE54002)", function() {
  banner("BRCA RNA EXTERNAL VALIDATION - GSE54002")
  es <- load_series_matrix("GSE54002")
  pd <- Biobase::pData(es)
  y  <- get_labels(pd, "tumor/nontumor:ch1",
                   c("Non-tumor" = "Normal", "Tumor" = "Early"), "GSE54002")
  X  <- collapse_to_panel_genes(es, brca_rna_feats)
  run_external("BRCA", "RNA", brca_rna_feats, brca_full, brca_train_patients,
               X, y, zscore = TRUE)
})

# ---- PRAD methylation: GSE76938
prad_meth_ext_result <- safe_run("PRAD methylation (GSE76938)", function() {
  banner("PRAD METHYLATION EXTERNAL VALIDATION - GSE76938")
  es <- load_series_matrix("GSE76938")
  pd <- Biobase::pData(es)
  y  <- get_labels(pd, "tissue:ch1",
                   c("prostate benign tissue" = "Normal", "prostate cancer tissue" = "Early"), "GSE76938")
  X  <- extract_meth_matrix(es, prad_meth_feats)
  run_external("PRAD", "methylation", prad_meth_feats, prad_full, prad_train_patients,
               X, y, zscore = FALSE)
})

# ---- PRAD RNA: GSE29079
prad_rna_ext_result <- safe_run("PRAD RNA (GSE29079)", function() {
  banner("PRAD RNA EXTERNAL VALIDATION - GSE29079")
  es <- load_series_matrix("GSE29079")
  pd <- Biobase::pData(es)
  y  <- get_labels(pd, "tissuer type:ch1",
                   c("Benign" = "Normal", "Tumor" = "Early"), "GSE29079")
  X  <- collapse_to_panel_genes(es, prad_rna_feats)
  run_external("PRAD", "RNA", prad_rna_feats, prad_full, prad_train_patients,
               X, y, zscore = TRUE)
})

# ---- 8.9 FINAL SUMMARY
internal_file <- file.path(cfg8e$step8_dir, "Step8_internal_validation_summary.rds")
internal_summary <- if (file.exists(internal_file)) readRDS(internal_file) else {
  warning("Internal validation summary not found: ", internal_file); NULL
}

get_internal <- function(cohort, what) {
  v <- tryCatch({
    s <- internal_summary[[cohort]]
    if (what == "cv") {
      cv <- s$cv_performance
      cv$mean[cv$metric == "AUC"]
    } else {
      s$holdout_test$auc
    }
  }, error = function(e) NA_real_)
  if (length(v) == 0) NA_real_ else as.numeric(v[1])
}

fmt_auc <- function(x) {
  a <- x$auc
  if (is.null(a) || is.na(a)) "NA (see error/warnings above)" else sprintf("%.3f", a)
}
fmt_num <- function(a) if (is.na(a)) "NA" else sprintf("%.3f", a)

final_validation_summary <- list(
  BRCA = list(
    internal = if (!is.null(internal_summary)) internal_summary$BRCA else NULL,
    external = list(methylation = brca_meth_ext_result, rna = brca_rna_ext_result,
                    protein = "EXCLUDED -- no viable independent external RPPA source found")
  ),
  PRAD = list(
    internal = if (!is.null(internal_summary)) internal_summary$PRAD else NULL,
    external = list(methylation = prad_meth_ext_result, rna = prad_rna_ext_result)
  )
)
saveRDS(final_validation_summary, file.path(cfg8e$out_dir, "Step8_FINAL_validation_summary.rds"))

summary_tbl <- data.frame(
  cohort = c("BRCA", "BRCA", "PRAD", "PRAD"),
  omics  = c("methylation", "RNA", "methylation", "RNA"),
  geo    = c("GSE60185", "GSE54002", "GSE76938", "GSE29079"),
  n      = c(brca_meth_ext_result$n, brca_rna_ext_result$n,
             prad_meth_ext_result$n, prad_rna_ext_result$n),
  auc    = c(brca_meth_ext_result$auc, brca_rna_ext_result$auc,
             prad_meth_ext_result$auc, prad_rna_ext_result$auc),
  features_used = c(
    paste(brca_meth_ext_result$features_used, collapse = ";"),
    paste(brca_rna_ext_result$features_used,  collapse = ";"),
    paste(prad_meth_ext_result$features_used, collapse = ";"),
    paste(prad_rna_ext_result$features_used,  collapse = ";")),
  error  = c(
    if (is.null(brca_meth_ext_result$error)) "" else brca_meth_ext_result$error,
    if (is.null(brca_rna_ext_result$error))  "" else brca_rna_ext_result$error,
    if (is.null(prad_meth_ext_result$error)) "" else prad_meth_ext_result$error,
    if (is.null(prad_rna_ext_result$error))  "" else prad_rna_ext_result$error),
  stringsAsFactors = FALSE
)
write.csv(summary_tbl, file.path(cfg8e$out_dir, "Step8_FINAL_external_summary.csv"), row.names = FALSE)

banner("STEP 8.9 FINAL VALIDATION SUMMARY")
cat("\nBRCA:\n")
cat("  Internal CV AUC (mean):", fmt_num(get_internal("BRCA", "cv")), "\n")
cat("  Internal held-out test AUC:", fmt_num(get_internal("BRCA", "holdout")), "\n")
cat("  External methylation AUC (GSE60185):", fmt_auc(brca_meth_ext_result), "\n")
cat("  External RNA AUC (GSE54002):", fmt_auc(brca_rna_ext_result), "\n")
cat("\nPRAD:\n")
cat("  Internal CV AUC (mean):", fmt_num(get_internal("PRAD", "cv")), "\n")
cat("  Internal held-out test AUC:", fmt_num(get_internal("PRAD", "holdout")), "\n")
cat("  External methylation AUC (GSE76938):", fmt_auc(prad_meth_ext_result), "\n")
cat("  External RNA AUC (GSE29079):", fmt_auc(prad_rna_ext_result), "\n")
cat("\nAll outputs written to:\n", cfg8e$out_dir, "\n")
cat("============================================================\n")
