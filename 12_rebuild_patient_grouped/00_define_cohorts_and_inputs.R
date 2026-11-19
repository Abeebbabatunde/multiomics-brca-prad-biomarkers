# 00_define_cohorts_and_inputs.R
# Step 12.0 (Steps A and B): defines the rebuild cohorts and builds ONE harmonised set of inputs for both cancers.
# It replaces 11/00_build_inputs.R for the rebuild and leaves every original file unchanged (outputs end in _rebuild).
#
# FROZEN DEFINITIONS (edit only before the first full run, then record any change):
#  Cohort, breast : early stage = AJCC stage I or II from the CLINICAL file (clinical_brca.csv, 804 patients), not from the RNA
#                   metadata (which lacks stage for 52 patients); tumour = sample code 01, normal = code 11; one sample per
#                   patient and class (first sorted barcode); the patient needs a tumour AND a normal sample in BOTH layers.
#  Cohort, prostate: unchanged from 05_0b (Gleason score <= 7; stage is not available); same pairing rules.
#  Methylation    : both cancers: the normalised beta matrices from Step 2; probes restricted to those present in BOTH cancers
#                   and complete in the cohort samples; beta clipped to 0.001-0.999, then M-values.
#  RNA            : both cancers: DESeq2 vst with blind = TRUE (labels NOT used) on the saved dds objects; genes with a raw count >= 10
#                   in >= 5 cohort samples and present in BOTH cancers.
# Protein is not part of the rebuild (30 breast patients only; the panels are rebuilt for methylation and RNA).
# Edit project_dir only.

project_dir <- "C:/Users/USER/Documents/MyThesis_MultiOmics"
setwd(project_dir)
suppressPackageStartupMessages({library(data.table); library(DESeq2)})
out_dir <- "12_rebuild_inputs"; dir.create(out_dir, showWarnings = FALSE)

fix   <- function(x) gsub(".", "-", x, fixed = TRUE)
code  <- function(b) substr(fix(b), 14, 15)
pat   <- function(b) substr(fix(b), 1, 12)
to_M  <- function(beta) { beta <- pmin(pmax(beta, 0.001), 0.999); log2(beta / (1 - beta)) }
pick  <- function(d, p, cd) sort(d$barcode[d$patient == p & d$code == cd])[1]
both_ok <- function(d, early) intersect(unique(d$patient[d$code == "11"]), unique(d$patient[d$code == "01" & d$patient %in% early]))

## ------------------------------ 1. cohorts ------------------------------
barcodes_of <- function(x) { x <- fix(x); x <- x[grepl("^TCGA", x)]; data.frame(barcode = x, patient = pat(x), code = code(x), stringsAsFactors = FALSE) }
meth_hdr <- function(f) barcodes_of(names(fread(f, nrows = 0)))

rm_b <- as.data.frame(readRDS("02_ProcessedData/Breast/brca_rna_meta.rds"))
rna_b <- barcodes_of(rm_b$barcode); meth_b <- meth_hdr("02_ProcessedData/Breast/brca_meth_normalized.csv")
clin_b <- read.csv("01_RawData/Breast/clinical_brca.csv", stringsAsFactors = FALSE)
early_b <- unique(clin_b$submitter_id[grepl("^Stage (I|II)[AB]?$", clin_b$ajcc_pathologic_stage)])
cat("Breast early-stage patients (clinical file):", length(early_b), "\n")
cat("Of these, also early stage in the RNA metadata:", sum(early_b %in% unique(rm_b$patient[grepl("^Stage (I|II)[AB]?$", rm_b$ajcc_pathologic_stage)])), "\n")

rm_p <- as.data.frame(readRDS("02_ProcessedData/Prostate/prad_rna_meta.rds"))
rna_p <- barcodes_of(rm_p$barcode); meth_p <- meth_hdr("02_ProcessedData/Prostate/prad_meth_normalized.csv")
early_p <- unique(rm_p$patient[!is.na(rm_p$gleason_score) & rm_p$gleason_score <= 7])

make_manifest <- function(rna, meth, early) {
  ps <- sort(intersect(both_ok(rna, early), both_ok(meth, early)))
  do.call(rbind, lapply(ps, function(p) do.call(rbind, lapply(c("Normal", "Early"), function(g) {
    cd <- ifelse(g == "Normal", "11", "01")
    data.frame(patient = p, group = g, rna = pick(rna, p, cd), meth = pick(meth, p, cd), stringsAsFactors = FALSE)
  }))))
}
man_b <- make_manifest(rna_b, meth_b, early_b); man_p <- make_manifest(rna_p, meth_p, early_p)
cat("Breast pairs (methylation + RNA):", length(unique(man_b$patient)), " [thesis: 51]\n")
cat("Prostate pairs (methylation + RNA):", length(unique(man_p$patient)), " [thesis: 17]\n")
write.csv(man_b, file.path(out_dir, "cohort_BRCA_rebuild.csv"), row.names = FALSE)
write.csv(man_p, file.path(out_dir, "cohort_PRAD_rebuild.csv"), row.names = FALSE)

## ------------------------------ 2. methylation (common probes) ------------------------------
read_beta <- function(f, needed) {                                  # reads only the cohort columns (the full matrices are large)
  hdr <- names(fread(f, nrows = 0)); id <- hdr[!grepl("^TCGA", fix(hdr))][1]
  sel <- c(id, hdr[fix(hdr) %in% fix(needed)])
  m <- fread(f, select = sel, data.table = FALSE)
  rn <- m[[1]]; m <- as.matrix(m[, -1, drop = FALSE]); colnames(m) <- fix(colnames(m)); rownames(m) <- rn; m
}
mb <- read_beta("02_ProcessedData/Breast/brca_meth_normalized.csv", man_b$meth)
mp <- read_beta("02_ProcessedData/Prostate/prad_meth_normalized.csv", man_p$meth)
common_probes <- intersect(rownames(mb), rownames(mp))
cat("Methylation probes: breast", nrow(mb), " prostate", nrow(mp), " common", length(common_probes), "\n")
build_meth <- function(m, man, tag) {
  stopifnot(all(fix(man$meth) %in% colnames(m)))
  mm <- m[common_probes, fix(man$meth)]; mm <- mm[rowSums(is.na(mm)) == 0, ]
  X <- t(to_M(mm))
  meta <- data.frame(sample = rownames(X), patient = man$patient, class = ifelse(man$group == "Normal", "Normal", "Tumour"))
  saveRDS(X, file.path(out_dir, paste0("X_methylation_", tag, "_rebuild.rds")))
  write.csv(meta, file.path(out_dir, paste0("meta_", tag, "_methylation_rebuild.csv")), row.names = FALSE)
  cat(tag, "methylation:", dim(X), "\n"); X
}
Xmb <- build_meth(mb, man_b, "BRCA"); Xmp <- build_meth(mp, man_p, "PRAD")
rm(mb, mp); gc()

## ------------------------------ 3. RNA (blind vst, common genes) ------------------------------
build_rna <- function(dds_file, man, tag) {
  dds <- readRDS(dds_file); colnames(dds) <- fix(colnames(dds))
  stopifnot(all(man$rna %in% colnames(dds)))
  raw <- counts(dds)[, man$rna, drop = FALSE]
  keep <- rowSums(raw >= 10) >= 5                                    # label-free rule applied to the cohort samples
  v <- assay(vst(dds, blind = TRUE))[, man$rna, drop = FALSE]        # blind = TRUE: tumour/normal labels are not used
  list(v = v, keep = rownames(dds)[keep])
}
rb <- build_rna("02_ProcessedData/Breast/dds_brca.rds", man_b, "BRCA")
rp <- build_rna("02_ProcessedData/Prostate/dds_prad.rds", man_p, "PRAD")
common_genes <- intersect(rb$keep, rp$keep)
cat("RNA genes kept: breast", length(rb$keep), " prostate", length(rp$keep), " common", length(common_genes), "\n")
save_rna <- function(v, man, tag) {
  X <- t(v[common_genes, , drop = FALSE]); X <- X[, apply(X, 2, function(z) !anyNA(z)), drop = FALSE]
  meta <- data.frame(sample = rownames(X), patient = man$patient, class = ifelse(man$group == "Normal", "Normal", "Tumour"))
  saveRDS(X, file.path(out_dir, paste0("X_rna_", tag, "_rebuild.rds")))
  write.csv(meta, file.path(out_dir, paste0("meta_", tag, "_rna_rebuild.csv")), row.names = FALSE)
  cat(tag, "RNA:", dim(X), "\n")
}
save_rna(rb$v, man_b, "BRCA"); save_rna(rp$v, man_p, "PRAD")
cat("\nDone. Check the pair counts against the thesis (51 and 17) and record any difference as a deliberate change of cohort definition.\n")
