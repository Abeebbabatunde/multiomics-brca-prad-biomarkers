# 04_biological_and_subtype_check.R
# Step 12.4 (Objective 2, step G): biological description and breast-subtype consistency of the REBUILT panels.
#  (1) Direction and size of the tumour-minus-normal difference of every rebuilt feature, paired within patient, in the rebuild cohort
#      (RNA: difference in variance-stabilised units, approximately log2 fold change; CpG: difference in beta value),
#      with the share of patients whose difference has the same sign as the median.
#  (2) CpG annotation (gene and region) from the Illumina 450K annotation package, as in 10_final_panel/10_3b.
#  (3) Breast only: the same paired difference within each PAM50 subtype of the tumour (subtype labels from brca_meth_clinical_full.rds, as in 09_1).
# The features were selected for this very contrast, so the p-values below are DESCRIPTIVE ONLY (circular) and are not evidence of discovery;
# the effect size, its direction and its consistency are what to read. Subtypes with few patients are flagged.
# Inputs: 00_define_cohorts_and_inputs.R and 01_stability_selection_rebuild.R. Prostate methylation has an empty rebuilt panel and is skipped. Edit project_dir only.

project_dir <- "C:/Users/USER/Documents/MyThesis_MultiOmics"
setwd(project_dir)

layers <- list(
  BRCA_methylation = list(x = "12_rebuild_inputs/X_methylation_BRCA_rebuild.rds", meta = "12_rebuild_inputs/meta_BRCA_methylation_rebuild.csv", type = "meth", cancer = "BRCA"),
  BRCA_RNA         = list(x = "12_rebuild_inputs/X_rna_BRCA_rebuild.rds",         meta = "12_rebuild_inputs/meta_BRCA_rna_rebuild.csv",         type = "rna",  cancer = "BRCA"),
  PRAD_methylation = list(x = "12_rebuild_inputs/X_methylation_PRAD_rebuild.rds", meta = "12_rebuild_inputs/meta_PRAD_methylation_rebuild.csv", type = "meth", cancer = "PRAD"),
  PRAD_RNA         = list(x = "12_rebuild_inputs/X_rna_PRAD_rebuild.rds",         meta = "12_rebuild_inputs/meta_PRAD_rna_rebuild.csv",         type = "rna",  cancer = "PRAD"))
min_n_subtype <- 8      # subtypes with fewer paired patients are flagged as low-n

## ---- breast PAM50 subtype of the tumour, by patient ----
fixb <- function(x) gsub(".", "-", x, fixed = TRUE)
mc <- as.data.frame(readRDS("02_ProcessedData/Breast/brca_meth_clinical_full.rds"))
bc <- fixb(mc$barcode)
pam <- data.frame(patient = substr(bc, 1, 12), type = substr(bc, 14, 15), pam50 = as.character(mc$paper_BRCA_Subtype_PAM50), stringsAsFactors = FALSE)
pam <- pam[pam$type == "01" & !is.na(pam$pam50) & nzchar(pam$pam50) & pam$pam50 != "NA", ]
pam <- pam[!duplicated(pam$patient), ]

## ---- paired tumour-minus-normal difference of one feature, one value per patient ----
paired_diff <- function(X, meta, f, type) {
  meta <- meta[match(rownames(X), meta$sample), ]
  v <- X[, f]; if (type == "meth") v <- 2^v / (1 + 2^v)               # M-value -> beta
  pts <- unique(meta$patient)
  d <- sapply(pts, function(p) v[which(meta$patient == p & meta$class == "Tumour")[1]] - v[which(meta$patient == p & meta$class == "Normal")[1]])
  setNames(as.numeric(d), pts)
}
summ <- function(d) {
  d <- d[!is.na(d)]
  data.frame(n_pairs = length(d), median_diff = median(d), share_same_direction = mean(sign(d) == sign(median(d))),
             p_descriptive = suppressWarnings(if (length(d) >= 3) wilcox.test(d)$p.value else NA_real_))
}

eff <- list(); sub <- list(); cpgs <- c()
for (nm in names(layers)) {
  cfg <- layers[[nm]]
  panel <- read.csv(paste0(nm, "_rebuilt_panel.csv"), stringsAsFactors = FALSE)
  if (nrow(panel) == 0) { cat(nm, ": empty rebuilt panel, skipped\n"); next }
  X <- readRDS(cfg$x); meta <- read.csv(cfg$meta, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(panel))) {
    f <- panel$feature[i]; lab <- if (!is.null(panel$symbol) && !is.na(panel$symbol[i])) panel$symbol[i] else f
    d <- paired_diff(X, meta, f, cfg$type)
    eff[[length(eff) + 1]] <- cbind(layer = nm, feature = f, label = lab, selection_frequency = panel$selection_frequency[i],
                                    measure = ifelse(cfg$type == "meth", "delta beta (tumour - normal)", "VST difference (tumour - normal)"), summ(d))
    if (cfg$cancer == "BRCA") {
      st <- pam$pam50[match(names(d), pam$patient)]
      for (s in sort(unique(na.omit(st)))) {
        r <- summ(d[which(st == s)])
        sub[[length(sub) + 1]] <- cbind(layer = nm, feature = f, label = lab, subtype = s, r, low_n = r$n_pairs < min_n_subtype)
      }
    }
    if (cfg$type == "meth") cpgs <- c(cpgs, f)
  }
}
eff <- do.call(rbind, eff); eff$q_descriptive <- p.adjust(eff$p_descriptive, "BH")
cat("\n(1) Paired tumour-minus-normal difference of every rebuilt feature (p-values descriptive only):\n")
print(eff[, c("layer", "label", "selection_frequency", "measure", "n_pairs", "median_diff", "share_same_direction", "p_descriptive")], row.names = FALSE, digits = 3)
write.csv(eff, "rebuild_04_effects.csv", row.names = FALSE)

sub <- do.call(rbind, sub)
cat("\n(3) Same difference within each PAM50 subtype (low_n = fewer than", min_n_subtype, "paired patients):\n")
print(sub[, c("layer", "label", "subtype", "n_pairs", "median_diff", "share_same_direction", "low_n")], row.names = FALSE, digits = 3)
write.csv(sub, "rebuild_04_breast_subtype.csv", row.names = FALSE)
cat("\nPAM50 subtype counts among the paired breast patients:\n"); print(table(pam$pam50[match(unique(read.csv("12_rebuild_inputs/meta_BRCA_rna_rebuild.csv")$patient), pam$patient)], useNA = "ifany"))

## ---- (2) CpG annotation ----
pk <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
if (length(cpgs) && requireNamespace(pk, quietly = TRUE)) {
  e <- new.env(); data("Other", package = pk, envir = e)
  an <- as.data.frame(e$Other[cpgs, c("UCSC_RefGene_Name", "UCSC_RefGene_Group")]); an$probe <- cpgs
  cat("\n(2) CpG annotation:\n"); print(an, row.names = FALSE)
  write.csv(an, "rebuild_04_cpg_annotation.csv", row.names = FALSE)
} else cat("\n(2) CpG annotation skipped: install", pk, "via BiocManager::install()\n")
