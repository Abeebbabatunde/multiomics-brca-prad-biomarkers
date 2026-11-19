# 05_external_validation_rebuilt_panels.R
# Step 12.5 (Objective 2, step H): external validation of the REBUILT panels, with the thesis external method unchanged
# (same GEO cohorts, same random forest, ntree = 1000, seed 42, per-dataset z-scoring for RNA, no external labels used for anything but scoring).
#   Breast methylation  : GSE60185 (rebuilt CpGs)          Breast RNA   : GSE54002 (rebuilt genes)
#   Prostate RNA        : GSE29079 (rebuilt genes)         Prostate methylation: no rebuilt panel, not tested.
# The model for each layer is trained on ALL patients of the rebuild cohort (the rebuilt panel was chosen on these patients, so the internal fit is
# not a test); the external cohort is used for scoring only. Features absent from an external platform are dropped and reported (disclose in Methods).
# Reports AUC with DeLong and bootstrap 95% CIs, and sensitivity / specificity at the model's own call (no external labels used to set a cut-off)
# with exact (Clopper-Pearson) intervals; the thesis results for the locked panels are printed beside them if the thesis file is found.
# Needs the series matrices already cached from the thesis run (08_FinalResults/External_GEO_cache) or an internet connection. Edit project_root only.

project_root <- "C:/Users/USER/Documents/MyThesis_MultiOmics"
setwd(project_root)
cfg8e <- list(geo_cache_dir = file.path(project_root, "08_FinalResults", "External_GEO_cache"),
              out_dir = file.path(project_root, "12_rebuild_external"), seed = 42, rf_ntree = 1000)
dir.create(cfg8e$out_dir, recursive = TRUE, showWarnings = FALSE); dir.create(cfg8e$geo_cache_dir, recursive = TRUE, showWarnings = FALSE)
options(timeout = max(600, getOption("timeout")))

## ---- helpers copied unchanged from 08_validation/08_6_external_validation.R (GEO download, probe-to-gene mapping, labels, z-score) ----
suppressPackageStartupMessages({
  library(Biobase)
  library(GEOquery)
  library(randomForest)
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


## ---- training data from the rebuild cohort, in the same units as the external data ----
panel_of <- function(nm) read.csv(paste0(nm, "_rebuilt_panel.csv"), stringsAsFactors = FALSE)
train_data <- function(nm, type) {
  tag <- sub("_.*", "", nm)
  if (type == "meth") {
    X <- readRDS(paste0("12_rebuild_inputs/X_methylation_", tag, "_rebuild.rds")); meta <- read.csv(paste0("12_rebuild_inputs/meta_", tag, "_methylation_rebuild.csv"), stringsAsFactors = FALSE)
    p <- panel_of(nm); X <- 2^X[, p$feature, drop = FALSE] / (1 + 2^X[, p$feature, drop = FALSE])      # rebuild matrices hold M-values; the GEO data hold betas
  } else {
    X <- readRDS(paste0("12_rebuild_inputs/X_rna_", tag, "_rebuild.rds")); meta <- read.csv(paste0("12_rebuild_inputs/meta_", tag, "_rna_rebuild.csv"), stringsAsFactors = FALSE)
    p <- panel_of(nm); p <- p[!is.na(p$symbol) & nzchar(p$symbol), ]
    X <- X[, p$feature, drop = FALSE]; colnames(X) <- p$symbol                                          # genes by symbol, as in the thesis external step
  }
  meta <- meta[match(rownames(X), meta$sample), ]
  list(X = as.data.frame(X), y = factor(ifelse(meta$class == "Normal", "Normal", "Early"), levels = c("Normal", "Early")))
}

## ---- train on the rebuild cohort, score the external cohort ----
ext_run <- function(nm, tr, ext_X, ext_y, zscore) {
  ext_y <- ext_y[rownames(ext_X)]; keep <- !is.na(ext_y); ext_X <- ext_X[keep, , drop = FALSE]; ext_y <- ext_y[keep]
  feats <- intersect(colnames(tr$X), colnames(ext_X)); miss <- setdiff(colnames(tr$X), feats)
  cat("\n", nm, ": features used", length(feats), "of", ncol(tr$X), if (length(miss)) paste0(" | not on the external platform: ", paste(miss, collapse = ", ")), "\n", sep = "")
  if (!length(feats)) stop("No rebuilt features available in the external data.")
  Xe <- as.data.frame(ext_X[, feats, drop = FALSE]); Xe[] <- lapply(Xe, as.numeric)
  ok <- complete.cases(Xe); Xe <- Xe[ok, , drop = FALSE]; ye <- factor(ext_y[ok], levels = c("Normal", "Early"))
  Xi <- tr$X[, feats, drop = FALSE]; Xi[] <- lapply(Xi, as.numeric)
  if (zscore) { Xi <- zscore_cols(Xi); Xe <- zscore_cols(Xe) }
  set.seed(cfg8e$seed)
  model <- randomForest::randomForest(x = Xi, y = tr$y, ntree = cfg8e$rf_ntree)
  prob <- predict(model, Xe, type = "prob")[, "Early"]; pred <- predict(model, Xe, type = "response")
  write.csv(data.frame(sample = rownames(Xe), true = as.character(ye), pred = as.character(pred), prob_early = as.numeric(prob)),
            file.path(cfg8e$out_dir, paste0(nm, "_rebuilt_external_predictions.csv")), row.names = FALSE)
  r <- pROC::roc(ye, prob, levels = c("Normal", "Early"), direction = "<", quiet = TRUE)
  cd <- as.numeric(pROC::ci.auc(r, method = "delong"))
  cb <- as.numeric(pROC::ci.auc(r, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE, progress = "none"))
  nt <- sum(ye == "Early"); nn <- sum(ye == "Normal"); tp <- sum(ye == "Early" & pred == "Early"); tn <- sum(ye == "Normal" & pred == "Normal")
  s <- binom.test(tp, nt)$conf.int; p <- binom.test(tn, nn)$conf.int
  data.frame(cohort = nm, n_features_used = length(feats), n_features_missing = length(miss), missing = paste(miss, collapse = ";"),
             n_tumour = nt, n_normal = nn, AUC = round(as.numeric(pROC::auc(r)), 3),
             DeLong_low = round(cd[1], 3), DeLong_high = round(cd[3], 3), boot_low = round(cb[1], 3), boot_high = round(cb[3], 3),
             sensitivity = round(tp / nt, 3), sens_low = round(s[1], 3), sens_high = round(s[2], 3),
             specificity = round(tn / nn, 3), spec_low = round(p[1], 3), spec_high = round(p[2], 3))
}
safe <- function(label, f) tryCatch(f(), error = function(e) { cat("\n!! ", label, " FAILED: ", conditionMessage(e), "\n", sep = ""); NULL })

res <- list()
res$BRCA_methylation <- safe("BRCA methylation", function() {
  es <- load_series_matrix("GSE60185"); pd <- Biobase::pData(es)
  y <- get_labels(pd, "disease state:ch1", c("normal" = "Normal", "breast cancer" = "Early"), "GSE60185")
  ext_run("BRCA_methylation", train_data("BRCA_methylation", "meth"), extract_meth_matrix(es, panel_of("BRCA_methylation")$feature), y, zscore = FALSE) })
res$BRCA_RNA <- safe("BRCA RNA", function() {
  es <- load_series_matrix("GSE54002"); pd <- Biobase::pData(es)
  y <- get_labels(pd, "tumor/nontumor:ch1", c("Non-tumor" = "Normal", "Tumor" = "Early"), "GSE54002")
  tr <- train_data("BRCA_RNA", "rna")
  ext_run("BRCA_RNA", tr, collapse_to_panel_genes(es, colnames(tr$X)), y, zscore = TRUE) })
res$BRCA_RNA_WISP1alias <- safe("BRCA RNA, CCN4 matched as WISP1 (sensitivity analysis)", function() {   # CCN4 is also known as WISP1; the older array may list it under that name
  es <- load_series_matrix("GSE54002"); pd <- Biobase::pData(es)
  y <- get_labels(pd, "tumor/nontumor:ch1", c("Non-tumor" = "Normal", "Tumor" = "Early"), "GSE54002")
  tr <- train_data("BRCA_RNA", "rna")
  ext <- collapse_to_panel_genes(es, sub("^CCN4$", "WISP1", colnames(tr$X)))
  colnames(ext)[colnames(ext) == "WISP1"] <- "CCN4"
  ext_run("BRCA_RNA_WISP1alias", tr, ext, y, zscore = TRUE) })
res$PRAD_RNA <- safe("PRAD RNA", function() {
  es <- load_series_matrix("GSE29079"); pd <- Biobase::pData(es)
  y <- get_labels(pd, "tissuer type:ch1", c("Benign" = "Normal", "Tumor" = "Early"), "GSE29079")
  tr <- train_data("PRAD_RNA", "rna")
  ext_run("PRAD_RNA", tr, collapse_to_panel_genes(es, colnames(tr$X)), y, zscore = TRUE) })
cat("\nPRAD_methylation: empty rebuilt panel, not tested.\n")

out <- do.call(rbind, res)
cat("\n=== Rebuilt panels, external validation ===\n"); print(out, row.names = FALSE)
write.csv(out, file.path(cfg8e$out_dir, "rebuilt_external_validation_with_CIs.csv"), row.names = FALSE)

lock <- file.path(project_root, "08_FinalResults", "Step8_External_Validation", "external_validation_with_CIs.csv")
if (file.exists(lock) && !is.null(out)) {
  lk <- read.csv(lock, stringsAsFactors = FALSE)
  cmp <- merge(out[, c("cohort", "n_features_used", "AUC", "DeLong_low", "DeLong_high")], lk[, c("cohort", "AUC", "DeLong_low", "DeLong_high")],
               by = "cohort", suffixes = c("_rebuilt", "_locked"))
  cat("\n=== Rebuilt vs locked panels (same external cohorts) ===\n"); print(cmp, row.names = FALSE)
  write.csv(cmp, file.path(cfg8e$out_dir, "rebuilt_vs_locked_external.csv"), row.names = FALSE)
} else cat("\nLocked-panel external results file not found; comparison skipped.\n")
