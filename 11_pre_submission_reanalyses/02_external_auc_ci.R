# 02_external_auc_ci.R
# FIX 2: 95% confidence intervals for the external AUCs, and exact intervals for sensitivity and specificity.
# Input: one prediction file per external cohort, written by the Step 8 external-validation scripts, with columns
#   true (Normal / Early), prob_early (probability from the LOCKED model), pred (call at the locked cut-off).
# "Early" is the tumour class. The Youden cut-off is calculated only to illustrate threshold transfer
# (it uses the external labels, so it is NOT a validation result).

suppressPackageStartupMessages(library(pROC))
set.seed(2026)
setwd("C:/Users/USER/Documents/MyThesis_MultiOmics")
d <- "08_FinalResults/Step8_External_Validation/"
cohorts <- c("BRCA_methylation", "BRCA_RNA", "PRAD_methylation", "PRAD_RNA")

res <- do.call(rbind, lapply(cohorts, function(nm) {
  x <- read.csv(paste0(d, nm, "_external_predictions.csv"), stringsAsFactors = FALSE)
  x$true <- factor(x$true, levels = c("Normal", "Early"))
  stopifnot(!anyNA(x$true))
  r    <- roc(x$true, x$prob_early, levels = c("Normal", "Early"), direction = "<", quiet = TRUE)
  ci_d <- as.numeric(ci.auc(r, method = "delong"))
  ci_b <- as.numeric(ci.auc(r, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE, progress = "none"))
  n_t <- sum(x$true == "Early");  n_n <- sum(x$true == "Normal")
  tp  <- sum(x$true == "Early"  & x$pred == "Early")
  tn  <- sum(x$true == "Normal" & x$pred == "Normal")
  s <- binom.test(tp, n_t)$conf.int                                   # Clopper-Pearson
  p <- binom.test(tn, n_n)$conf.int
  yo <- coords(r, "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)[1, ]
  data.frame(cohort = nm, n_tumour = n_t, n_normal = n_n,
             AUC = round(as.numeric(auc(r)), 3),
             DeLong_low = round(ci_d[1], 3), DeLong_high = round(ci_d[3], 3),
             boot_low = round(ci_b[1], 3), boot_high = round(ci_b[3], 3),
             sensitivity = round(tp / n_t, 3), sens_low = round(s[1], 3), sens_high = round(s[2], 3),
             specificity = round(tn / n_n, 3), spec_low = round(p[1], 3), spec_high = round(p[2], 3),
             posthoc_threshold = round(yo$threshold, 3),
             posthoc_sens = round(yo$sensitivity, 3), posthoc_spec = round(yo$specificity, 3))
}))
print(res, row.names = FALSE)
write.csv(res, paste0(d, "external_validation_with_CIs.csv"), row.names = FALSE)
