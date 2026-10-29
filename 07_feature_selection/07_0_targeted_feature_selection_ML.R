# ============================================================
# STEP 7.0 — TARGETED BIOMARKER FEATURE SELECTION & MACHINE LEARNING
# (patient-level split, bootstrap stability selection, patient-grouped CV)
# Assembled from three pasted segments: config + packages, helpers + BRCA build,
# PRAD build + selection function + run. Code unchanged.
# ============================================================

# ------------------------------------------------------------
# 0. CONFIG — edit paths/params here, nothing else in the script
# ------------------------------------------------------------
project_root <- "C:/Users/USER/Documents/MyThesis_MultiOmics"

cfg <- list(
  # Step 6 consensus gene panels (targeted feature universe)
  brca_consensus_file = file.path(project_root, "06_FinalResults", "Step6_Network",
                                   "BRCA_3omics", "consensus_genes.csv"),
  prad_consensus_file = file.path(project_root, "06_FinalResults", "Step6_Network",
                                   "PRAD_meth_rna", "consensus_genes.csv"),

  # BRCA processed matrices + label manifest + RPPA protein source
  brca_meth_file  = file.path(project_root, "02_ProcessedData", "Breast", "brca_meth_normalized.csv"),
  brca_rna_file   = file.path(project_root, "02_ProcessedData", "Breast", "brca_rna_normalized.csv"),
  brca_label_file = file.path(project_root, "00.5_clinical_metadata_harmonization",
                               "sample_manifest_paired_normal_early.csv"),
  brca_rppa_pairing_file = file.path(project_root, "01_RawData", "Breast",
                                      "brca_protein_sample_pairing.csv"),
  brca_rppa_dir = file.path(project_root, "GDCdata", "TCGA-BRCA", "Proteome_Profiling",
                             "Protein_Expression_Quantification"),
  # Authoritative 16-patient 3-omics Normal/Early manifest (patient, group,
  # rna, meth, prot columns all in one file -- this IS the official
  # 3omics_16pairs cohort used by the Step 5 DIABLO/MOFA model)
  brca_3omics_label_file = file.path(project_root, "05_FinalResults", "Step5.0_Integration",
                                      "BRCA", "cohort_3omics.csv"),

  # PRAD processed matrices + OFFICIAL 17-pair label manifest
  prad_meth_file  = file.path(project_root, "02_ProcessedData", "Prostate", "prad_meth_normalized.csv"),
  prad_rna_file   = file.path(project_root, "02_ProcessedData", "Prostate", "prad_rna_normalized.csv"),
  prad_label_file = file.path(project_root, "05_FinalResults", "Step5.0_Integration",
                               "PRAD", "cohort_meth_rna.csv"),

  # Output location for this step
  out_dir = file.path(project_root, "07_FinalResults", "Step7_FeatureSelection_ML"),

  # ML params
  seed              = 42,
  test_fraction     = 0.25,   # held-out test set fraction (final, untouched evaluation)
  cv_folds          = 5,      # cross-validation folds for tuning/RFE
  rf_ntree          = 1000,
  lasso_alpha       = 1,      # 1 = LASSO, 0 = ridge, between = elastic net
  min_method_votes  = 2,      # a feature must be picked by >= this many of
                               # {LASSO, Boruta/RF, SVM-RFE} to enter final panel

  # Stability selection params (7.4) — repeat selection across bootstraps of
  # the TRAINING data only; a feature must be picked frequently, not just once
  n_bootstrap         = 50,   # number of bootstrap resamples of the train set
  stability_threshold = 0.6,  # a feature must be selected in >= this fraction
                               # of bootstraps (per method) to count as "stable"
  boruta_maxRuns      = 50,   # reduced from Boruta's default 100 for speed,
                               # since this runs inside every bootstrap

  # Repeated CV performance assessment (7.6) — instead of trusting one
  # small holdout, get a distribution of performance across many resamples
  cv_repeats          = 10    # repeats of cv_folds-fold CV for final panel evaluation
)

dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------
required_pkgs <- c("data.table", "dplyr", "glmnet", "Boruta", "randomForest",
                    "caret", "e1071", "org.Hs.eg.db", "AnnotationDbi")
to_install <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) {
  cat("Installing missing packages:", paste(to_install, collapse = ", "), "\n")
  cran_pkgs <- setdiff(to_install, c("org.Hs.eg.db", "AnnotationDbi"))
  if (length(cran_pkgs) > 0) install.packages(cran_pkgs)
  bioc_pkgs <- intersect(to_install, c("org.Hs.eg.db", "AnnotationDbi"))
  if (length(bioc_pkgs) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(bioc_pkgs, update = FALSE, ask = FALSE)
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(glmnet)
  library(Boruta)
  library(randomForest)
  library(caret)
  library(e1071)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(dplyr)  # reload last so dplyr::select/filter/mutate/rename aren't
                  # masked by AnnotationDbi/org.Hs.eg.db's own S4 generics
})

set.seed(cfg$seed)

# ------------------------------------------------------------
# 2. HELPER FUNCTIONS
# ------------------------------------------------------------

# TCGA barcode sample-type code: positions 14-15 ("01"=tumor,"11"=normal, etc.)
sample_type_code <- function(barcode) substr(barcode, 14, 15)
# TCGA patient ID: first 3 barcode fields
patient_id_from_barcode <- function(barcode) substr(barcode, 1, 12)

# Read a features-as-rows omics CSV (first column = feature ID) as a matrix,
# subset to the requested feature IDs, return samples-as-rows matrix.
load_and_subset_omics <- function(path, feature_ids, feature_col_name = "V1") {
  dt <- fread(path, check.names = FALSE)
  setnames(dt, 1, "feature_id")
  dt_sub <- dt[feature_id %in% feature_ids]
  if (nrow(dt_sub) == 0) {
    warning("No matching features found in ", path)
  }
  feat_ids <- dt_sub$feature_id
  mat <- as.matrix(dt_sub[, -1, with = FALSE])
  rownames(mat) <- feat_ids
  t(mat)  # return samples x features
}

# Map RNA gene symbols (from consensus panel) to Ensembl gene IDs (no version),
# so they can be matched to the versioned Ensembl row names in the RNA matrix.
symbols_to_ensembl <- function(symbols) {
  map <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = unique(symbols), keytype = "SYMBOL", columns = "ENSEMBL"
  ))
  map <- map[!is.na(map$ENSEMBL), ]
  map
}

# Strip Ensembl version suffix (".15" etc.) from matrix column/row names
strip_ensembl_version <- function(x) sub("\\..*$", "", x)

# ------------------------------------------------------------
# 3. BUILD BRCA TARGETED SAMPLE x FEATURE MATRIX (METH + RNA + RPPA PROTEIN)
# ------------------------------------------------------------
cat("\n=== Building BRCA targeted matrix (true 3-omics, 16-patient cohort) ===\n")

brca_consensus <- read.csv(cfg$brca_consensus_file, stringsAsFactors = FALSE)
brca_meth_genes  <- brca_consensus %>% filter(block == "meth") %>% pull(feature) %>% unique()
brca_rna_symbols <- brca_consensus %>% filter(block == "rna")  %>% pull(gene)    %>% unique()
brca_prot_feats  <- brca_consensus %>% filter(block == "prot") %>% pull(feature) %>% unique()

brca_labels_all <- read.csv(cfg$brca_label_file, stringsAsFactors = FALSE)
stopifnot(all(c("patient", "group", "rna", "meth") %in% names(brca_labels_all)))

# --- Use the AUTHORITATIVE 16-patient 3-omics manifest directly ---
# cohort_3omics.csv already has patient, group, rna, meth, prot (short barcode)
# all in one file -- no need to reconstruct the cohort from separate sources.
brca_labels <- read.csv(cfg$brca_3omics_label_file, stringsAsFactors = FALSE)
stopifnot(all(c("patient", "group", "rna", "meth", "prot") %in% names(brca_labels)))
cat("BRCA 3-omics official manifest:", n_distinct(brca_labels$patient), "patients,",
    nrow(brca_labels), "sample rows. Class balance:\n")
print(table(brca_labels$group))

# --- Methylation: subset directly by CpG probe ID, restrict to cohort ---
brca_meth_mat <- load_and_subset_omics(cfg$brca_meth_file, brca_meth_genes)
brca_meth_mat <- brca_meth_mat[rownames(brca_meth_mat) %in% brca_labels$meth, , drop = FALSE]
cat("BRCA meth matrix (targeted, 3-omics cohort):", dim(brca_meth_mat)[1], "samples x",
    dim(brca_meth_mat)[2], "features (of", length(brca_meth_genes), "requested)\n")

# --- RNA: map symbols -> Ensembl, subset by stripped Ensembl ID, restrict to cohort ---
sym_map <- symbols_to_ensembl(brca_rna_symbols)
cat("Mapped", length(unique(sym_map$SYMBOL)), "of", length(brca_rna_symbols),
    "BRCA consensus RNA symbols to Ensembl IDs\n")

brca_rna_dt <- fread(cfg$brca_rna_file, check.names = FALSE)
setnames(brca_rna_dt, 1, "feature_id")
brca_rna_dt[, ensembl_noversion := strip_ensembl_version(feature_id)]
brca_rna_dt_sub <- brca_rna_dt[ensembl_noversion %in% sym_map$ENSEMBL]
brca_rna_dt_sub <- merge(brca_rna_dt_sub,
                          unique(sym_map[, c("ENSEMBL", "SYMBOL")]),
                          by.x = "ensembl_noversion", by.y = "ENSEMBL", all.x = TRUE)

rna_feat_names <- brca_rna_dt_sub$SYMBOL
rna_cols <- setdiff(names(brca_rna_dt_sub), c("feature_id", "ensembl_noversion", "SYMBOL"))
brca_rna_mat <- t(as.matrix(brca_rna_dt_sub[, ..rna_cols]))
colnames(brca_rna_mat) <- rna_feat_names
brca_rna_mat <- brca_rna_mat[rownames(brca_rna_mat) %in% brca_labels$rna, , drop = FALSE]
cat("BRCA rna matrix (targeted, 3-omics cohort):", dim(brca_rna_mat)[1], "samples x",
    dim(brca_rna_mat)[2], "features (of", length(brca_rna_symbols), "requested)\n")

# --- Protein: read RPPA files for the 16-patient cohort, subset to consensus features ---
# brca_labels$prot holds the exact short barcode (e.g. "TCGA-BH-A0HA-11A") used
# to identify each sample's RPPA file. RPPA file naming on disk is the FULL
# barcode (e.g. TCGA-BH-A0HA-11A-31D-A12R-05_RPPA_data.tsv), so we match on
# the short-barcode prefix.
all_rppa_files <- list.files(cfg$brca_rppa_dir, recursive = TRUE,
                              pattern = "RPPA_data\\.tsv$", full.names = TRUE)
rppa_barcode_from_path <- sub("_RPPA_data\\.tsv$", "", basename(all_rppa_files))
rppa_lookup <- data.frame(file = all_rppa_files, full_barcode = rppa_barcode_from_path,
                           short_barcode = substr(rppa_barcode_from_path, 1, 16),
                           stringsAsFactors = FALSE)

read_rppa_targeted <- function(file_path, keep_targets) {
  dt <- fread(file_path, select = c("peptide_target", "protein_expression"))
  dt <- dt[peptide_target %in% keep_targets]
  setNames(as.list(dt$protein_expression), dt$peptide_target)
}

prot_rows <- list()
for (i in seq_len(nrow(brca_labels))) {
  short_bc <- brca_labels$prot[i]
  match_row <- rppa_lookup[rppa_lookup$short_barcode == short_bc, ]
  if (nrow(match_row) == 0) {
    warning("No RPPA file found for ", short_bc, " -- skipping")
    next
  }
  vals <- read_rppa_targeted(match_row$file[1], brca_prot_feats)
  vals$prot <- short_bc
  prot_rows[[length(prot_rows) + 1]] <- vals
}
brca_prot_df <- bind_rows(lapply(prot_rows, as.data.frame, stringsAsFactors = FALSE))

cat("BRCA protein matrix (targeted, RPPA):", nrow(brca_prot_df), "samples x",
    length(brca_prot_feats), "features (of", length(brca_prot_feats), "requested)\n")

saveRDS(brca_prot_df, file.path(cfg$out_dir, "BRCA_RPPA_protein_targeted.rds"))

# --- Merge meth + rna + protein into one sample x feature table via patient+group ---
brca_meth_df <- data.frame(meth = rownames(brca_meth_mat), brca_meth_mat,
                            check.names = FALSE, row.names = NULL)
brca_rna_df  <- data.frame(rna = rownames(brca_rna_mat), brca_rna_mat,
                            check.names = FALSE, row.names = NULL)

brca_full <- brca_labels %>%
  dplyr::select(patient, group, rna, meth, prot) %>%
  inner_join(brca_meth_df, by = "meth") %>%
  inner_join(brca_rna_df,  by = "rna") %>%
  inner_join(brca_prot_df, by = "prot")

cat("BRCA final targeted 3-omics dataset:", nrow(brca_full), "samples x",
    ncol(brca_full) - 4, "features (meth+rna+prot). Class balance:\n")
print(table(brca_full$group))

saveRDS(brca_full, file.path(cfg$out_dir, "BRCA_targeted_dataset.rds"))
write.csv(brca_full, file.path(cfg$out_dir, "BRCA_targeted_dataset.csv"), row.names = FALSE)

# ------------------------------------------------------------
# 4. BUILD PRAD TARGETED SAMPLE x FEATURE MATRIX (METH + RNA)
# ------------------------------------------------------------
cat("\n=== Building PRAD targeted matrix (official 17-pair manifest) ===\n")

prad_consensus <- read.csv(cfg$prad_consensus_file, stringsAsFactors = FALSE)
prad_meth_probes <- prad_consensus %>% filter(block == "meth") %>% pull(feature) %>% unique()
prad_rna_symbols <- prad_consensus %>% filter(block == "rna") %>% pull(gene) %>% unique()

# Official, pipeline-validated PRAD Normal/Early manifest (17 patients x 2 groups),
# same structure as the BRCA manifest: patient, group, rna, meth columns with
# exact TCGA sample barcodes. This is NOT hand-built -- it is the same manifest
# your Step 5 DIABLO/MOFA integration for PRAD meth_rna_17pairs used.
prad_labels <- read.csv(cfg$prad_label_file, stringsAsFactors = FALSE)
stopifnot(all(c("patient", "group", "rna", "meth") %in% names(prad_labels)))
cat("PRAD official manifest:", n_distinct(prad_labels$patient), "patients,",
    nrow(prad_labels), "sample rows. Class balance:\n")
print(table(prad_labels$group))

# --- Methylation subset ---
prad_meth_mat <- load_and_subset_omics(cfg$prad_meth_file, prad_meth_probes)
prad_meth_mat <- prad_meth_mat[rownames(prad_meth_mat) %in% prad_labels$meth, , drop = FALSE]
cat("PRAD meth matrix (targeted):", dim(prad_meth_mat)[1], "samples x",
    dim(prad_meth_mat)[2], "features (of", length(prad_meth_probes), "requested)\n")

# --- RNA subset (symbol -> Ensembl mapping again) ---
prad_sym_map <- symbols_to_ensembl(prad_rna_symbols)
prad_rna_dt <- fread(cfg$prad_rna_file, check.names = FALSE)
setnames(prad_rna_dt, 1, "feature_id")
prad_rna_dt[, ensembl_noversion := strip_ensembl_version(feature_id)]
prad_rna_dt_sub <- prad_rna_dt[ensembl_noversion %in% prad_sym_map$ENSEMBL]
prad_rna_dt_sub <- merge(prad_rna_dt_sub,
                          unique(prad_sym_map[, c("ENSEMBL", "SYMBOL")]),
                          by.x = "ensembl_noversion", by.y = "ENSEMBL", all.x = TRUE)

prad_rna_feat_names <- prad_rna_dt_sub$SYMBOL
prad_rna_cols_data <- setdiff(names(prad_rna_dt_sub), c("feature_id", "ensembl_noversion", "SYMBOL"))
prad_rna_mat <- t(as.matrix(prad_rna_dt_sub[, ..prad_rna_cols_data]))
colnames(prad_rna_mat) <- prad_rna_feat_names
prad_rna_mat <- prad_rna_mat[rownames(prad_rna_mat) %in% prad_labels$rna, , drop = FALSE]
cat("PRAD rna matrix (targeted):", dim(prad_rna_mat)[1], "samples x",
    dim(prad_rna_mat)[2], "features (of", length(prad_rna_symbols), "requested)\n")

# --- Merge PRAD meth + rna ---
prad_meth_df <- data.frame(meth = rownames(prad_meth_mat), prad_meth_mat,
                            check.names = FALSE, row.names = NULL)
prad_rna_df  <- data.frame(rna = rownames(prad_rna_mat), prad_rna_mat,
                            check.names = FALSE, row.names = NULL)

prad_full <- prad_labels %>%
  dplyr::select(patient, group, rna, meth) %>%
  inner_join(prad_meth_df, by = "meth") %>%
  inner_join(prad_rna_df,  by = "rna")

cat("PRAD final targeted dataset:", nrow(prad_full), "samples x",
    ncol(prad_full) - 4, "features. Class balance:\n")
print(table(prad_full$group))

saveRDS(prad_full, file.path(cfg$out_dir, "PRAD_targeted_dataset.rds"))
write.csv(prad_full, file.path(cfg$out_dir, "PRAD_targeted_dataset.csv"), row.names = FALSE)

# ------------------------------------------------------------
# 5. LEAKAGE-SAFE, STABILITY-SELECTED FEATURE SELECTION + ML
#    (reusable for BRCA and PRAD)
#
# 7.2 Leakage-safe filtering: near-zero-variance filter is fit on the TRAIN
#     split only (never on test), then the same retained feature set is
#     applied to test. Nothing about the test set influences any filtering
#     or selection decision anywhere in this function.
# PATIENT-LEVEL INTEGRITY (applies throughout 7.2-7.6): this dataset has 2
#     rows per patient (paired Normal+Early samples). Every split in this
#     function -- the initial train/test split, every bootstrap resample in
#     stability selection, and every CV fold in performance assessment -- is
#     done by PATIENT, not by row, so a patient's paired samples always stay
#     together on one side. A row-level split here would let the model be
#     evaluated on a patient it effectively already trained on.
# 7.3 Feature selection: LASSO, Boruta(RF), SVM-RFE -- as before, but now
#     run inside a bootstrap loop rather than once (see 7.4).
# 7.4 Stable features: the whole 7.3 selection process is repeated across
#     cfg$n_bootstrap patient-level bootstrap resamples of the TRAINING data.
#     A feature only counts as selected by a method if it was picked in >=
#     cfg$stability_threshold fraction of the bootstraps for that method.
#     The final consensus panel requires >= cfg$min_method_votes of the
#     3 methods' STABLE feature sets to agree, not just one lucky split.
# 7.5 ML classification: Random Forest only (by design decision -- with
#     n=24-26 training samples, comparing 3 classifiers on one holdout adds
#     noise without adding real evidence; one well cross-validated model is
#     more defensible here than three under-powered ones).
# 7.6 Performance assessment: repeated, patient-grouped k-fold CV (cfg$cv_repeats
#     x cfg$cv_folds, via caret::groupKFold) on the training data for a
#     distribution of AUC/accuracy, PLUS a single untouched final holdout-test
#     evaluation as a secondary, independent check (not used to tune anything).
# 7.7 Preliminary biomarker signature: the final stable consensus panel.
#
# KNOWN LIMITATION (disclose in thesis, not fixable without redoing Steps 5-6):
#     the Step 6 consensus gene panel itself was built using DIABLO/MOFA on
#     the FULL cohort (all 16 BRCA / 17 PRAD patients) before any train/test
#     split existed here. So the candidate feature LIST was chosen with
#     knowledge of every patient, including whoever ends up in this step's
#     test set. Everything inside this function (selection, stability,
#     performance) is done train-only and patient-safe, but the upstream
#     panel definition itself is not nested -- a fully leakage-free design
#     would require re-running Steps 5-6 inside every CV fold, which is a
#     much larger undertaking than this step. This is a common, accepted
#     limitation of two-stage discovery-then-validate designs in small-cohort
#     biomarker studies, and should be named as such in the discussion.
# ------------------------------------------------------------
run_feature_selection_and_ml <- function(full_df, id_cols, label_col = "group",
                                          cohort_name = "COHORT", cfg) {

  cat("\n============================================================\n")
  cat("FEATURE SELECTION + ML —", cohort_name, "\n")
  cat("============================================================\n")

  feature_cols <- setdiff(names(full_df), id_cols)
  X_all <- full_df[, feature_cols, drop = FALSE]
  X_all[] <- lapply(X_all, as.numeric)
  y_all <- factor(full_df[[label_col]], levels = c("Normal", "Early"))

  # --- PATIENT-LEVEL split (not sample-level) ---
  # Each patient contributes 2 rows (Normal + Early paired samples). Splitting
  # by row would let one patient's Normal land in train while their Early
  # lands in test (or vice versa) -- the model would effectively be evaluated
  # on someone it already trained on. So we split by unique patient ID first,
  # then take BOTH of that patient's rows into whichever side they land on.
  stopifnot("patient" %in% names(full_df))
  patients <- unique(full_df$patient)
  set.seed(cfg$seed)
  n_test_patients <- max(1, round(length(patients) * cfg$test_fraction))
  test_patients  <- sample(patients, n_test_patients)
  train_patients <- setdiff(patients, test_patients)

  train_mask <- full_df$patient %in% train_patients
  test_mask  <- full_df$patient %in% test_patients

  X_train_raw <- X_all[train_mask, , drop = FALSE]; y_train <- y_all[train_mask]
  X_test_raw  <- X_all[test_mask, , drop = FALSE];  y_test  <- y_all[test_mask]
  train_patient_ids <- full_df$patient[train_mask]  # needed for patient-level bootstrap below

  cat("Patient-level split:", length(train_patients), "train patients (",
      sum(train_mask), "samples ) |", length(test_patients), "test patients (",
      sum(test_mask), "samples ). No patient appears in both.\n")

  # --- 7.2 Leakage-safe near-zero-variance filter: fit on TRAIN only ---
  nzv <- caret::nearZeroVar(X_train_raw)
  if (length(nzv) > 0) {
    cat("Dropping", length(nzv), "near-zero-variance features (fit on train only)\n")
    X_train <- X_train_raw[, -nzv, drop = FALSE]
    X_test  <- X_test_raw[, -nzv, drop = FALSE]
  } else {
    X_train <- X_train_raw; X_test <- X_test_raw
  }

  cat("Train:", nrow(X_train), "| Test:", nrow(X_test),
      "| Features entering selection:", ncol(X_train), "\n")

  # ---------------- 7.3 + 7.4: bootstrap-repeated feature selection ----------------
  cat("\n--- Stability selection across", cfg$n_bootstrap, "bootstraps (patient-level) ---\n")

  # Patient-level bootstrap: resample unique TRAIN patients with replacement
  # (never touches test patients), then pull in BOTH rows for every sampled
  # patient. This keeps each patient's paired Normal+Early samples together
  # in every bootstrap draw, consistent with the paired study design.
  unique_train_patients <- unique(train_patient_ids)
  set.seed(cfg$seed)
  boot_idx_list <- lapply(seq_len(cfg$n_bootstrap), function(i) {
    sampled_patients <- sample(unique_train_patients, length(unique_train_patients), replace = TRUE)
    unlist(lapply(sampled_patients, function(p) which(train_patient_ids == p)))
  })

  all_feature_names <- colnames(X_train)
  lasso_hits <- setNames(rep(0L, length(all_feature_names)), all_feature_names)
  rf_hits    <- setNames(rep(0L, length(all_feature_names)), all_feature_names)
  svm_hits   <- setNames(rep(0L, length(all_feature_names)), all_feature_names)

  for (b in seq_along(boot_idx_list)) {
    idx <- boot_idx_list[[b]]
    Xb <- X_train[idx, , drop = FALSE]
    yb <- y_train[idx]
    if (length(unique(yb)) < 2) next  # skip degenerate resamples

    # LASSO on this bootstrap
    lasso_b_feats <- tryCatch({
      cvb <- glmnet::cv.glmnet(as.matrix(Xb), yb, family = "binomial",
                                alpha = cfg$lasso_alpha, nfolds = min(cfg$cv_folds, 3))
      cb <- coef(cvb, s = "lambda.1se")
      setdiff(rownames(cb)[which(cb[, 1] != 0)], "(Intercept)")
    }, error = function(e) character(0))
    lasso_hits[lasso_b_feats] <- lasso_hits[lasso_b_feats] + 1L

    # Boruta on this bootstrap (fewer maxRuns for speed inside the loop)
    rf_b_feats <- tryCatch({
      bo <- Boruta::Boruta(x = Xb, y = yb, ntree = cfg$rf_ntree,
                            maxRuns = cfg$boruta_maxRuns, doTrace = 0)
      bo <- Boruta::TentativeRoughFix(bo)
      names(bo$finalDecision)[bo$finalDecision == "Confirmed"]
    }, error = function(e) character(0))
    rf_hits[rf_b_feats] <- rf_hits[rf_b_feats] + 1L

    # SVM-RFE on this bootstrap
    svm_b_feats <- tryCatch({
      sizes_b <- unique(round(seq(2, ncol(Xb), length.out = min(6, ncol(Xb)))))
      rfe_ctrl_b <- caret::rfeControl(functions = caret::caretFuncs, method = "cv",
                                       number = min(cfg$cv_folds, 3), verbose = FALSE)
      rfe_b <- caret::rfe(Xb, yb, sizes = sizes_b, rfeControl = rfe_ctrl_b, method = "svmLinear")
      caret::predictors(rfe_b)
    }, error = function(e) character(0))
    svm_hits[svm_b_feats] <- svm_hits[svm_b_feats] + 1L

    if (b %% 10 == 0) cat("  bootstrap", b, "/", cfg$n_bootstrap, "done\n")
  }

  freq_df <- data.frame(
    feature      = all_feature_names,
    lasso_freq   = lasso_hits / cfg$n_bootstrap,
    rf_freq      = rf_hits    / cfg$n_bootstrap,
    svm_freq     = svm_hits   / cfg$n_bootstrap,
    stringsAsFactors = FALSE
  )

  stable_lasso <- freq_df$feature[freq_df$lasso_freq >= cfg$stability_threshold]
  stable_rf    <- freq_df$feature[freq_df$rf_freq    >= cfg$stability_threshold]
  stable_svm   <- freq_df$feature[freq_df$svm_freq   >= cfg$stability_threshold]

  cat("\nStable features (selected in >=", cfg$stability_threshold * 100, "% of bootstraps):\n")
  cat("  LASSO:", length(stable_lasso), "| Boruta/RF:", length(stable_rf),
      "| SVM-RFE:", length(stable_svm), "\n")

  # ---------------- Final consensus panel from STABLE sets ----------------
  all_feats <- union(union(stable_lasso, stable_rf), stable_svm)
  vote_df <- freq_df[freq_df$feature %in% all_feats, ]
  vote_df$in_lasso <- vote_df$feature %in% stable_lasso
  vote_df$in_rf    <- vote_df$feature %in% stable_rf
  vote_df$in_svm   <- vote_df$feature %in% stable_svm
  vote_df$n_votes  <- rowSums(vote_df[, c("in_lasso", "in_rf", "in_svm")])
  vote_df <- vote_df[order(-vote_df$n_votes, -vote_df$lasso_freq - vote_df$rf_freq - vote_df$svm_freq), ]
  consensus_features <- vote_df$feature[vote_df$n_votes >= cfg$min_method_votes]

  cat("\n7.7 Preliminary biomarker signature (>=", cfg$min_method_votes,
      "of 3 methods, stability-filtered):", length(consensus_features), "features\n")
  print(consensus_features)

  write.csv(freq_df, file.path(cfg$out_dir, paste0(cohort_name, "_bootstrap_selection_frequencies.csv")),
            row.names = FALSE)
  write.csv(vote_df, file.path(cfg$out_dir, paste0(cohort_name, "_feature_selection_votes.csv")),
            row.names = FALSE)

  # ---------------- 7.6 Performance assessment ----------------
  cv_summary <- NULL
  cm <- NULL
  auc_test <- NA

  if (length(consensus_features) >= 2) {
    final_train <- X_train[, consensus_features, drop = FALSE]
    final_test  <- X_test[, consensus_features, drop = FALSE]

    # (a) Repeated k-fold CV on the training data -- distribution of performance.
    # Folds are grouped by PATIENT (via caret::groupKFold), repeated cfg$cv_repeats
    # times with different random fold assignments, so a patient's paired
    # Normal+Early rows are never split across train/validation within a fold.
    cv_index <- list()
    for (r in seq_len(cfg$cv_repeats)) {
      set.seed(cfg$seed + r)
      folds_r <- caret::groupKFold(train_patient_ids, k = cfg$cv_folds)
      names(folds_r) <- paste0("Repeat", r, ".", names(folds_r))
      cv_index <- c(cv_index, folds_r)
    }
    ctrl_rep <- caret::trainControl(method = "cv", index = cv_index,
                                     classProbs = TRUE,
                                     summaryFunction = caret::twoClassSummary,
                                     savePredictions = "final")
    set.seed(cfg$seed)
    rf_model <- caret::train(x = final_train, y = y_train, method = "rf",
                              trControl = ctrl_rep, metric = "ROC")

    cv_summary <- rf_model$resample  # one row per fold x repeat
    cat("\nRepeated CV performance (", cfg$cv_repeats, "repeats x", cfg$cv_folds, "folds ):\n")
    cat("  Mean ROC (AUC):", round(mean(cv_summary$ROC), 3),
        "  SD:", round(sd(cv_summary$ROC), 3), "\n")
    cat("  Mean Sens:", round(mean(cv_summary$Sens), 3),
        "  Mean Spec:", round(mean(cv_summary$Spec), 3), "\n")

    # (b) Independent final holdout test -- untouched throughout, secondary check
    pred <- predict(rf_model, final_test)
    pred_prob <- predict(rf_model, final_test, type = "prob")
    cm <- caret::confusionMatrix(pred, y_test, positive = "Early")
    roc_obj <- pROC::roc(response = y_test, predictor = pred_prob$Early, quiet = TRUE)
    auc_test <- as.numeric(pROC::auc(roc_obj))

    cat("\nFinal held-out test performance (secondary check, RF on consensus panel):\n")
    print(cm)
    cat("Held-out test AUC:", auc_test, "\n")

    saveRDS(rf_model, file.path(cfg$out_dir, paste0(cohort_name, "_final_RF_model.rds")))
    write.csv(cv_summary, file.path(cfg$out_dir, paste0(cohort_name, "_repeated_cv_performance.csv")),
              row.names = FALSE)
  } else {
    cat("\nFewer than 2 consensus features — skipping classifier evaluation.\n")
  }

  list(stable_lasso = stable_lasso, stable_rf = stable_rf, stable_svm = stable_svm,
       consensus = consensus_features, votes = vote_df, freq = freq_df,
       cv_summary = cv_summary, confusion_matrix = cm, holdout_auc = auc_test)
}

if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
library(pROC)

brca_results <- run_feature_selection_and_ml(
  brca_full, id_cols = c("patient", "group", "rna", "meth", "prot"),
  cohort_name = "BRCA", cfg = cfg
)

prad_results <- run_feature_selection_and_ml(
  prad_full, id_cols = c("patient", "group", "rna", "meth"),
  cohort_name = "PRAD", cfg = cfg
)

# ------------------------------------------------------------
# 6. SAVE STEP 7 SUMMARY
# ------------------------------------------------------------
saveRDS(list(BRCA = brca_results, PRAD = prad_results),
        file.path(cfg$out_dir, "Step7_feature_selection_results.rds"))

cat("\n============================================================\n")
cat("STEP 7.0 COMPLETE. Outputs written to:\n", cfg$out_dir, "\n")
cat("============================================================\n")
