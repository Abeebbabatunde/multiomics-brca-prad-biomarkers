## 05_1a  BRCA integration: data preparation, MOFA2 and DIABLO (Step 5.0)
## Inputs : cohort_3omics.csv and cohort_meth_rna.csv (from 05_0a_BRCA_cohort.R)
## Outputs: 05_FinalResults/Step5.0_Integration/BRCA/{meth_rna_51pairs,3omics_16pairs}/
## Python : MOFA2 was run with use_basilisk = FALSE and mofapy2 0.7.5 under Python 3.11.9
##          (see 05_6_mofa_python_environment.R). Set RETICULATE_PYTHON before loading MOFA2/reticulate.
## Note   : an earlier version defined a wrapper called run_mofa, which masked MOFA2::run_mofa
##          ("unused arguments"); it was renamed do_mofa (earlier versions are in archive/).
Sys.setenv(RETICULATE_PYTHON = "C:/Users/USER/AppData/Local/Programs/Python/Python311/python.exe")
setwd("C:/Users/USER/Documents/MyThesis_MultiOmics")
suppressPackageStartupMessages({ library(data.table); library(DESeq2); library(MOFA2)
  library(mixOmics); library(org.Hs.eg.db); library(ggplot2) })
set.seed(123)
base <- "05_FinalResults/Step5.0_Integration/BRCA"
fix  <- function(x) gsub(".", "-", x, fixed = TRUE)
top_var <- function(m, n) m[order(apply(m, 1, var), decreasing = TRUE)[1:min(n, nrow(m))], , drop = FALSE]

c3 <- fread(file.path(base, "cohort_3omics.csv"))
c2 <- fread(file.path(base, "cohort_meth_rna.csv"))

## ---- RNA (VST) -> gene symbols ----
R   <- SummarizedExperiment::assay(readRDS("02_ProcessedData/Breast/brca_vst.rds"))[, unique(c2$rna)]
sym <- mapIds(org.Hs.eg.db, sub("\\..*", "", rownames(R)), "SYMBOL", "ENSEMBL", multiVals = "first")
R <- R[!is.na(sym), ]; sym <- sym[!is.na(sym)]
o <- order(apply(R, 1, var), decreasing = TRUE); R <- R[o, ]; sym <- sym[o]
R <- R[!duplicated(sym), ]; rownames(R) <- sym[!duplicated(sym)]

## ---- Methylation: read only cohort columns, beta -> M-values ----
f   <- "02_ProcessedData/Breast/brca_meth_normalized.csv"
hdr <- names(fread(f, nrows = 0)); idx <- match(unique(c2$meth), fix(hdr)); stopifnot(!anyNA(idx))
M   <- fread(f, select = c(1, idx)); pr <- M[[1]]
M   <- as.matrix(M[, -1, with = FALSE]); rownames(M) <- pr; colnames(M) <- fix(colnames(M))
M   <- M[complete.cases(M), unique(c2$meth)]
M   <- pmin(pmax(M, 1e-3), 1 - 1e-3); Mv <- log2(M / (1 - M)); rm(M); gc()
Mv  <- top_var(Mv, 50000)

## ---- Protein: long -> wide ----
pl <- fread("02_Normalized/Breast/brca_protein_normalized.csv"); pl[, bc := fix(sample_barcode)]
Pm <- dcast(pl, peptide_target ~ bc, value.var = "value")
rn <- Pm$peptide_target; Pm <- as.matrix(Pm[, -1, with = FALSE]); rownames(Pm) <- rn

## ---- assemble a cohort: samples named <patient>_N / <patient>_T ----
assemble <- function(coh, use_prot) {
  coh <- as.data.frame(coh)
  sid <- paste0(coh$patient, "_", ifelse(coh$group == "Normal", "N", "T"))
  B <- list(meth = Mv[, coh$meth], rna = R[, coh$rna])
  if (use_prot) B$prot <- Pm[, coh$prot]
  B <- lapply(B, function(m) { colnames(m) <- sid; m })
  list(X = B, group = factor(coh$group, levels = c("Normal", "Early")), patient = coh$patient, sid = sid)
}
ntop <- c(meth = 5000, rna = 5000, prot = 457)

## ---- MOFA ----
do_mofa <- function(A, tag, nf = 6) {
  od <- file.path(base, tag); dir.create(od, showWarnings = FALSE)
  X  <- mapply(top_var, A$X, ntop[names(A$X)], SIMPLIFY = FALSE)
  mo <- create_mofa(lapply(X, as.matrix))
  dopt <- get_default_data_options(mo); dopt$scale_views <- TRUE
  mopt <- get_default_model_options(mo); mopt$num_factors <- nf
  topt <- get_default_training_options(mo); topt$seed <- 123; topt$convergence_mode <- "slow"
  mo <- prepare_mofa(mo, data_options = dopt, model_options = mopt, training_options = topt)
  mo <- MOFA2::run_mofa(mo, outfile = file.path(od, "MOFA_model.hdf5"), use_basilisk = FALSE)
  samples_metadata(mo) <- data.frame(sample = A$sid, condition = A$group, patient = A$patient)
  ggsave(file.path(od, "MOFA_variance.png"), plot_variance_explained(mo), width = 6, height = 4)
  ggsave(file.path(od, "MOFA_factors.png"),
         plot_factor(mo, factors = 1:3, color_by = "condition", add_violin = TRUE), width = 7, height = 5)
  Z <- get_factors(mo)[[1]]; pts <- unique(A$patient)
  pv <- apply(Z, 2, function(z) wilcox.test(z[paste0(pts, "_T")], z[paste0(pts, "_N")], paired = TRUE)$p.value)
  cat("\nMOFA", tag, "- paired Wilcoxon p per factor (Early vs Normal):\n"); print(signif(pv, 3))
  for (v in names(X)) write.csv(get_weights(mo, views = v, as.data.frame = TRUE),
                                file.path(od, sprintf("MOFA_weights_%s.csv", v)), row.names = FALSE)
  invisible(mo)
}

## ---- DIABLO (paired: within-patient adjustment) ----
run_diablo <- function(A, tag, folds, grid) {
  od <- file.path(base, tag); dir.create(od, showWarnings = FALSE)
  X <- mapply(top_var, A$X, ntop[names(A$X)], SIMPLIFY = FALSE)
  X <- lapply(X, function(m) withinVariation(t(m), design = data.frame(sample = A$patient)))
  Y <- A$group
  design <- matrix(0.1, length(X), length(X), dimnames = list(names(X), names(X))); diag(design) <- 0
  tn <- tune.block.splsda(X, Y, ncomp = 2, test.keepX = grid[names(X)], design = design,
                          validation = "Mfold", folds = folds, nrepeat = 5, dist = "centroids.dist")
  dia <- block.splsda(X, Y, ncomp = 2, keepX = tn$choice.keepX, design = design)
  pf  <- perf(dia, validation = "Mfold", folds = folds, nrepeat = 10, dist = "centroids.dist")
  cat("\nDIABLO", tag, "- weighted-vote error (optimistic, small n):\n"); print(pf$WeightedVote.error.rate)
  png(file.path(od, "DIABLO_indiv.png"),    2000, 1600, res = 250); plotIndiv(dia, ellipse = TRUE, legend = TRUE); dev.off()
  png(file.path(od, "DIABLO_diablo.png"),   2000, 1600, res = 250); plotDiablo(dia, ncomp = 1); dev.off()
  try({ png(file.path(od, "DIABLO_circos.png"), 2000, 2000, res = 250); circosPlot(dia, cutoff = 0.6, line = TRUE); dev.off() })
  try({ png(file.path(od, "DIABLO_loadings.png"), 2400, 1800, res = 250); plotLoadings(dia, comp = 1, contrib = "max", method = "median"); dev.off() })
  sel <- do.call(rbind, lapply(names(X), function(b) do.call(rbind, lapply(1:2, function(k) {
    sv <- selectVar(dia, block = b, comp = k)[[b]]
    data.frame(block = b, comp = k, feature = sv$name, loading = sv$value$value.var) }))))
  write.csv(sel, file.path(od, "DIABLO_selected_features.csv"), row.names = FALSE)
  dmp <- fread("03_Results/Epigenomics/brca_DMPs.csv")
  cat("Selected methylation probes that are also DMPs:",
      round(mean(sel$feature[sel$block == "meth"] %in% dmp$Probe), 2), "\n")
  saveRDS(list(model = dia, tune = tn, perf = pf), file.path(od, "DIABLO_objects.rds"))
  invisible(dia)
}

## ---- RUN ----
A2 <- assemble(c2, FALSE)   # meth + RNA, 51 pairs (main)
A3 <- assemble(c3, TRUE)    # meth + RNA + protein, 16 pairs (exploratory)
g2 <- list(meth = c(10, 25, 50), rna = c(10, 25, 50))
g3 <- list(meth = c(10, 25, 50), rna = c(10, 25, 50), prot = c(5, 10, 20))

mofa2 <- do_mofa(A2, "meth_rna_51pairs", nf = 8)
dia2  <- run_diablo(A2, "meth_rna_51pairs", folds = 5, grid = g2)
mofa3 <- do_mofa(A3, "3omics_16pairs",  nf = 5)
dia3  <- run_diablo(A3, "3omics_16pairs",  folds = 4, grid = g3)

save(A2, A3, g2, g3, ntop, top_var, base, run_diablo, do_mofa,
     file = file.path(base, "prep_objects.RData"))
