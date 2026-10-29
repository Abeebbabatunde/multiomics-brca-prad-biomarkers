# ============================================================
# STEP 8.0-8.4 — LOCK THE MODEL + INTERNAL VALIDATION
# Rebuilt to follow this exact workflow:
#
#   8.0 Lock the model      -- freeze feature panel + algorithm/hyperparams
#                               + reproduce Step 7's exact patient split
#   8.1 Cross-validation /
#       resampling           -- repeated patient-grouped k-fold CV using
#                               ONLY the locked features, refit per fold
#                               (never re-selecting features)
#   8.2 Stability of
#       biomarkers            -- bootstrap the locked panel's feature
#                               importance (do locked features stay
#                               influential across resamples?)
#   8.3 Performance
#       estimation            -- summarize 8.1's resampled performance
#   8.4 Locked held-out
#       test evaluation       -- ONE final fixed model, trained on all
#                               train patients, scored ONCE on Step 7's
#                               untouched test patients
#
# Nothing here re-runs LASSO/Boruta/SVM-RFE. The feature panel is DONE as
# of Step 7 -- this script only asks "how well does that exact panel,
# with a fixed algorithm, generalize?"
# ============================================================

# ------------------------------------------------------------
# 0. CONFIG
# ------------------------------------------------------------
project_root <- "C:/Users/USER/Documents/MyThesis_MultiOmics"

cfg8 <- list(
  step7_dir = file.path(project_root, "07_FinalResults", "Step7_FeatureSelection_ML"),
  out_dir   = file.path(project_root, "08_FinalResults", "Step8_Validation"),

  seed          = 42,           # MUST match Step 7's seed to reproduce the exact same split
  test_fraction = 0.25,         # MUST match Step 7's test_fraction
  rf_ntree      = 1000,         # fixed algorithm hyperparameter -- not tuned per fold
  cv_folds      = 5,
  cv_repeats    = 10,
  n_bootstrap   = 50            # for 8.2 biomarker stability only
)

dir.create(cfg8$out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(caret)
  library(randomForest)
  library(pROC)
})

# ------------------------------------------------------------
# 8.0 LOCK THE MODEL
# ------------------------------------------------------------
lock_the_model <- function(full_df, id_cols, locked_panel, label_col = "group",
                            cohort_name, cfg) {

  cat("\n============================================================\n")
  cat("8.0 LOCK THE MODEL —", cohort_name, "\n")
  cat("============================================================\n")

  feature_cols <- setdiff(names(full_df), id_cols)
  stopifnot(all(locked_panel %in% feature_cols))

  X_all <- full_df[, locked_panel, drop = FALSE]
  X_all[] <- lapply(X_all, as.numeric)
  y_all <- factor(full_df[[label_col]], levels = c("Normal", "Early"))

  # Reproduce Step 7's EXACT patient-level split (same seed, same fraction,
  # same patient order) -- this is not a new split, it is the same one.
  patients <- unique(full_df$patient)
  set.seed(cfg$seed)
  n_test_patients <- max(1, round(length(patients) * cfg$test_fraction))
  test_patients  <- sample(patients, n_test_patients)
  train_patients <- setdiff(patients, test_patients)

  train_mask <- full_df$patient %in% train_patients
  test_mask  <- full_df$patient %in% test_patients

  cat("Locked feature panel (", length(locked_panel), "features ):",
      paste(locked_panel, collapse = ", "), "\n")
  cat("Locked algorithm: Random Forest, ntree =", cfg$rf_ntree, "(fixed, not tuned per fold)\n")
  cat("Reproduced Step 7 split:", length(train_patients), "train patients (",
      sum(train_mask), "samples ) |", length(test_patients), "test patients (",
      sum(test_mask), "samples )\n")
  cat("VERIFY this matches Step 7's printed split counts for", cohort_name,
      "before trusting downstream results.\n")

  write.csv(data.frame(patient = train_patients, split = "train"),
            file.path(cfg$out_dir, paste0(cohort_name, "_locked_train_patients.csv")), row.names = FALSE)
  write.csv(data.frame(patient = test_patients, split = "test"),
            file.path(cfg$out_dir, paste0(cohort_name, "_locked_test_patients.csv")), row.names = FALSE)

  list(
    cohort_name = cohort_name,
    panel = locked_panel,
    X_all = X_all, y_all = y_all,
    full_df = full_df,
    train_mask = train_mask, test_mask = test_mask,
    train_patients = train_patients, test_patients = test_patients
  )
}

# ------------------------------------------------------------
# 8.1 CROSS-VALIDATION / RESAMPLING (locked features, refit per fold)
# ------------------------------------------------------------
run_locked_cv <- function(locked, cfg) {
  cat("\n--- 8.1 Cross-validation / resampling —", locked$cohort_name,
      "(", cfg$cv_repeats, "repeats x", cfg$cv_folds, "folds ) ---\n")

  X_train <- locked$X_all[locked$train_mask, , drop = FALSE]
  y_train <- locked$y_all[locked$train_mask]
  train_patient_ids <- locked$full_df$patient[locked$train_mask]

  repeat_metrics <- list()

  for (r in seq_len(cfg$cv_repeats)) {
    set.seed(cfg$seed + r)
    k <- min(cfg$cv_folds, length(unique(train_patient_ids)))
    folds <- caret::groupKFold(train_patient_ids, k = k)  # list of TRAIN indices per fold

    oof_pred <- rep(NA_character_, length(y_train))
    oof_prob <- rep(NA_real_, length(y_train))

    for (f in seq_along(folds)) {
      tr_idx <- folds[[f]]
      te_idx <- setdiff(seq_len(nrow(X_train)), tr_idx)
      if (length(te_idx) == 0 || length(unique(y_train[tr_idx])) < 2) next

      set.seed(cfg$seed + r * 100 + f)
      rf_fold <- randomForest::randomForest(x = X_train[tr_idx, , drop = FALSE],
                                             y = y_train[tr_idx], ntree = cfg$rf_ntree)
      prob <- predict(rf_fold, X_train[te_idx, , drop = FALSE], type = "prob")[, "Early"]
      pred <- predict(rf_fold, X_train[te_idx, , drop = FALSE])
      oof_prob[te_idx] <- prob
      oof_pred[te_idx] <- as.character(pred)
    }

    valid <- !is.na(oof_pred)
    true_r <- y_train[valid]
    pred_r <- factor(oof_pred[valid], levels = c("Normal", "Early"))
    prob_r <- oof_prob[valid]

    cm_r <- caret::confusionMatrix(pred_r, true_r, positive = "Early")
    auc_r <- tryCatch(as.numeric(pROC::auc(pROC::roc(true_r, prob_r, quiet = TRUE))), error = function(e) NA)

    repeat_metrics[[r]] <- data.frame(
      repeat_num = r, auc = auc_r,
      accuracy = as.numeric(cm_r$overall["Accuracy"]),
      sensitivity = as.numeric(cm_r$byClass["Sensitivity"]),
      specificity = as.numeric(cm_r$byClass["Specificity"])
    )
  }

  metrics_df <- bind_rows(repeat_metrics)
  cat("Per-repeat pooled out-of-fold performance (each repeat covers every train sample once):\n")
  print(metrics_df)

  write.csv(metrics_df, file.path(cfg8$out_dir, paste0(locked$cohort_name, "_8.1_CV_repeat_metrics.csv")),
            row.names = FALSE)

  metrics_df
}

# ------------------------------------------------------------
# 8.2 STABILITY OF BIOMARKERS (locked panel's importance across bootstraps)
# ------------------------------------------------------------
run_biomarker_stability <- function(locked, cfg) {
  cat("\n--- 8.2 Stability of biomarkers —", locked$cohort_name,
      "(", cfg$n_bootstrap, "patient-level bootstraps ) ---\n")

  X_train <- locked$X_all[locked$train_mask, , drop = FALSE]
  y_train <- locked$y_all[locked$train_mask]
  train_patient_ids <- locked$full_df$patient[locked$train_mask]
  unique_train_patients <- unique(train_patient_ids)

  imp_matrix <- matrix(NA_real_, nrow = cfg$n_bootstrap, ncol = ncol(X_train),
                        dimnames = list(NULL, colnames(X_train)))

  set.seed(cfg$seed)
  for (b in seq_len(cfg$n_bootstrap)) {
    sampled_patients <- sample(unique_train_patients, length(unique_train_patients), replace = TRUE)
    idx <- unlist(lapply(sampled_patients, function(p) which(train_patient_ids == p)))
    Xb <- X_train[idx, , drop = FALSE]; yb <- y_train[idx]
    if (length(unique(yb)) < 2) next

    rf_b <- tryCatch(
      randomForest::randomForest(x = Xb, y = yb, ntree = cfg$rf_ntree, importance = TRUE),
      error = function(e) NULL)
    if (is.null(rf_b)) next
    imp_matrix[b, ] <- randomForest::importance(rf_b, type = 2)[, 1]  # MeanDecreaseGini
  }

  mean_imp <- colMeans(imp_matrix, na.rm = TRUE)
  sd_imp   <- apply(imp_matrix, 2, sd, na.rm = TRUE)
  # top-half ranking frequency: in what fraction of bootstraps does each
  # feature rank in the top half of the panel by importance?
  top_half_freq <- sapply(seq_len(ncol(imp_matrix)), function(j) {
    ranks <- apply(imp_matrix, 1, function(row) rank(-row, na.last = "keep")[j])
    mean(ranks <= ceiling(ncol(imp_matrix) / 2), na.rm = TRUE)
  })

  stability_df <- data.frame(
    feature = colnames(X_train),
    mean_importance = mean_imp,
    sd_importance = sd_imp,
    top_half_frequency = top_half_freq
  )
  stability_df <- stability_df[order(-stability_df$mean_importance), ]

  cat("Locked-panel feature importance stability across bootstraps:\n")
  print(stability_df)

  write.csv(stability_df, file.path(cfg8$out_dir, paste0(locked$cohort_name, "_8.2_biomarker_stability.csv")),
            row.names = FALSE)

  stability_df
}

# ------------------------------------------------------------
# 8.3 PERFORMANCE ESTIMATION (summarize 8.1)
# ------------------------------------------------------------
summarize_performance <- function(cv_metrics_df, cohort_name) {
  cat("\n--- 8.3 Performance estimation —", cohort_name, "---\n")

  summary_df <- data.frame(
    metric = c("AUC", "Accuracy", "Sensitivity", "Specificity"),
    mean = c(mean(cv_metrics_df$auc, na.rm = TRUE), mean(cv_metrics_df$accuracy, na.rm = TRUE),
             mean(cv_metrics_df$sensitivity, na.rm = TRUE), mean(cv_metrics_df$specificity, na.rm = TRUE)),
    sd = c(sd(cv_metrics_df$auc, na.rm = TRUE), sd(cv_metrics_df$accuracy, na.rm = TRUE),
           sd(cv_metrics_df$sensitivity, na.rm = TRUE), sd(cv_metrics_df$specificity, na.rm = TRUE))
  )
  cat("Internal cross-validated performance estimate (across", nrow(cv_metrics_df), "repeats):\n")
  print(summary_df)

  write.csv(summary_df, file.path(cfg8$out_dir, paste0(cohort_name, "_8.3_performance_estimate.csv")),
            row.names = FALSE)

  summary_df
}

# ------------------------------------------------------------
# 8.4 LOCKED HELD-OUT TEST EVALUATION
#     ONE fixed model, trained once on ALL train patients, scored ONCE.
# ------------------------------------------------------------
run_locked_holdout_test <- function(locked, cfg) {
  cat("\n--- 8.4 Locked held-out test evaluation —", locked$cohort_name, "---\n")

  X_train <- locked$X_all[locked$train_mask, , drop = FALSE]
  y_train <- locked$y_all[locked$train_mask]
  X_test  <- locked$X_all[locked$test_mask, , drop = FALSE]
  y_test  <- locked$y_all[locked$test_mask]

  set.seed(cfg$seed)
  final_locked_model <- randomForest::randomForest(x = X_train, y = y_train, ntree = cfg$rf_ntree)

  pred <- predict(final_locked_model, X_test)
  prob <- predict(final_locked_model, X_test, type = "prob")[, "Early"]

  cm <- caret::confusionMatrix(pred, y_test, positive = "Early")
  auc_val <- tryCatch(as.numeric(pROC::auc(pROC::roc(y_test, prob, quiet = TRUE))), error = function(e) NA)

  cat("Final locked-model performance on the never-touched held-out test set (",
      length(y_test), "samples ):\n")
  print(cm)
  cat("Held-out test AUC:", auc_val, "\n")

  saveRDS(final_locked_model, file.path(cfg8$out_dir, paste0(locked$cohort_name, "_8.4_LOCKED_model.rds")))
  write.csv(data.frame(true = as.character(y_test), pred = as.character(pred), prob_early = prob),
            file.path(cfg8$out_dir, paste0(locked$cohort_name, "_8.4_holdout_predictions.csv")),
            row.names = FALSE)

  list(model = final_locked_model, confusion_matrix = cm, auc = auc_val)
}

# ------------------------------------------------------------
# RUN 8.0-8.4 FOR BOTH COHORTS
# ------------------------------------------------------------
step7_results <- readRDS(file.path(cfg8$step7_dir, "Step7_feature_selection_results.rds"))
brca_full <- readRDS(file.path(cfg8$step7_dir, "BRCA_targeted_dataset.rds"))
prad_full <- readRDS(file.path(cfg8$step7_dir, "PRAD_targeted_dataset.rds"))

brca_locked <- lock_the_model(brca_full, id_cols = c("patient", "group", "rna", "meth", "prot"),
                               locked_panel = step7_results$BRCA$consensus,
                               cohort_name = "BRCA", cfg = cfg8)
prad_locked <- lock_the_model(prad_full, id_cols = c("patient", "group", "rna", "meth"),
                               locked_panel = step7_results$PRAD$consensus,
                               cohort_name = "PRAD", cfg = cfg8)

brca_cv       <- run_locked_cv(brca_locked, cfg8)
brca_stab     <- run_biomarker_stability(brca_locked, cfg8)
brca_perf     <- summarize_performance(brca_cv, "BRCA")
brca_holdout  <- run_locked_holdout_test(brca_locked, cfg8)

prad_cv       <- run_locked_cv(prad_locked, cfg8)
prad_stab     <- run_biomarker_stability(prad_locked, cfg8)
prad_perf     <- summarize_performance(prad_cv, "PRAD")
prad_holdout  <- run_locked_holdout_test(prad_locked, cfg8)

# ------------------------------------------------------------
# 8.9 (partial) — INTERNAL VALIDATION SUMMARY
# External (8.5-8.8) results get merged into this same object once the
# external validation script (Step8.5_External_Validation_*.R) is run.
# ------------------------------------------------------------
internal_validation_summary <- list(
  BRCA = list(cv_performance = brca_perf, biomarker_stability = brca_stab,
              holdout_test = list(confusion_matrix = brca_holdout$confusion_matrix,
                                   auc = brca_holdout$auc)),
  PRAD = list(cv_performance = prad_perf, biomarker_stability = prad_stab,
              holdout_test = list(confusion_matrix = prad_holdout$confusion_matrix,
                                   auc = prad_holdout$auc)),
  external = NULL  # filled in by the external validation script
)

saveRDS(internal_validation_summary, file.path(cfg8$out_dir, "Step8_internal_validation_summary.rds"))

cat("\n============================================================\n")
cat("STEP 8.0-8.4 COMPLETE. Locked models + all outputs written to:\n", cfg8$out_dir, "\n")
cat("Next: run the external validation inspection script, then apply the\n")
cat("LOCKED models saved here (", "BRCA_8.4_LOCKED_model.rds / PRAD_8.4_LOCKED_model.rds", ") to\n")
cat("the external cohort for steps 8.5-8.8.\n")
cat("============================================================\n")
