# STEP 9.2 SURVIVAL (BRCA overall survival) - PART 1/2
# Run AFTER Parts 1 and 2 of the 9.1 script, same R session
suppressPackageStartupMessages(library(survival))
if (!exists("build_dataset") || !exists("P")) stop("Run Parts 1 and 2 of the 9.1 script first.")
out_dir <- P("09_FinalResults", "Step9_Clinical")
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
feats <- ds$features[ds$features$omics != "protein", ]   # VAV1 (n=76) excluded
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
