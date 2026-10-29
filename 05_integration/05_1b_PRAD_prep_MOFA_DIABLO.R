## 05_1b  PRAD integration: data preparation, MOFA2 and DIABLO (Step 5.0)
## Inputs: cohort_meth_rna.csv (05_0b_PRAD_cohort.R). Python: see 05_6_mofa_python_environment.R
Sys.setenv(RETICULATE_PYTHON = "C:/Users/USER/AppData/Local/Programs/Python/Python311/python.exe")
setwd("C:/Users/USER/Documents/MyThesis_MultiOmics")
suppressPackageStartupMessages({ library(data.table); library(MOFA2); library(mixOmics); library(org.Hs.eg.db); library(ggplot2) })
set.seed(123)
base <- "05_FinalResults/Step5.0_Integration/PRAD"
fix  <- function(x) gsub(".", "-", x, fixed = TRUE)
top_var <- function(m, n) m[order(apply(m, 1, var), decreasing = TRUE)[1:min(n, nrow(m))], , drop = FALSE]

coh <- fread(file.path(base, "cohort_meth_rna.csv"))

## ---- RNA: read only cohort columns, map to gene symbols ----
f   <- "02_ProcessedData/Prostate/prad_rna_normalized.csv"
hdr <- names(fread(f, nrows = 0)); idx <- match(unique(coh$rna), fix(hdr)); stopifnot(!anyNA(idx))
Rd  <- fread(f, select = c(1, idx)); rn <- Rd[[1]]
R   <- as.matrix(Rd[, -1, with = FALSE]); rownames(R) <- rn; colnames(R) <- fix(colnames(R))
sym <- mapIds(org.Hs.eg.db, sub("\\..*", "", rownames(R)), "SYMBOL", "ENSEMBL", multiVals = "first")
R <- R[!is.na(sym), ]; sym <- sym[!is.na(sym)]
o <- order(apply(R, 1, var), decreasing = TRUE); R <- R[o, , drop = FALSE]; sym <- sym[o]
R <- R[!duplicated(sym), , drop = FALSE]; rownames(R) <- sym[!duplicated(sym)]

## ---- Methylation: read only cohort columns, beta -> M-values ----
f2  <- "02_ProcessedData/Prostate/prad_meth_normalized.csv"
hdr2 <- names(fread(f2, nrows = 0)); idx2 <- match(unique(coh$meth), fix(hdr2)); stopifnot(!anyNA(idx2))
Md  <- fread(f2, select = c(1, idx2)); pr <- Md[[1]]
M   <- as.matrix(Md[, -1, with = FALSE]); rownames(M) <- pr; colnames(M) <- fix(colnames(M))
M   <- M[complete.cases(M), unique(coh$meth), drop = FALSE]
M   <- pmin(pmax(M, 1e-3), 1 - 1e-3); Mv <- log2(M / (1 - M)); rm(M); gc()

## ---- assemble ----
sid <- paste0(coh$patient, "_", ifelse(coh$group == "Normal", "N", "T"))
A <- list(X = list(meth = Mv[, coh$meth], rna = R[, coh$rna]),
          group = factor(coh$group, levels = c("Normal", "Early")),
          patient = coh$patient, sid = sid)
A$X <- lapply(A$X, function(m) { colnames(m) <- sid; m })
ntop <- c(meth = 5000, rna = 5000)

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

## ---- DIABLO ----
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
  cat("\nDIABLO", tag, "- weighted-vote error:\n"); print(pf$WeightedVote.error.rate)
  png(file.path(od, "DIABLO_indiv.png"),  2000, 1600, res = 250); plotIndiv(dia, ellipse = TRUE, legend = TRUE); dev.off()
  png(file.path(od, "DIABLO_diablo.png"), 2000, 1600, res = 250); plotDiablo(dia, ncomp = 1); dev.off()
  try({ png(file.path(od, "DIABLO_circos.png"), 2000, 2000, res = 250); circosPlot(dia, cutoff = 0.6, line = TRUE); dev.off() })
  try({ png(file.path(od, "DIABLO_loadings.png"), 2400, 1800, res = 250); plotLoadings(dia, comp = 1, contrib = "max", method = "median"); dev.off() })
  sel <- do.call(rbind, lapply(names(X), function(b) do.call(rbind, lapply(1:2, function(k) {
    sv <- selectVar(dia, block = b, comp = k)[[b]]
    data.frame(block = b, comp = k, feature = sv$name, loading = sv$value$value.var) }))))
  write.csv(sel, file.path(od, "DIABLO_selected_features.csv"), row.names = FALSE)
  saveRDS(list(model = dia, tune = tn, perf = pf), file.path(od, "DIABLO_objects.rds"))
  invisible(dia)
}

g2 <- list(meth = c(10, 25, 50), rna = c(10, 25, 50))
save(A, ntop, top_var, base, run_diablo, do_mofa, file = file.path(base, "prep_objects.RData"))

dia  <- run_diablo(A, "meth_rna_17pairs", folds = 4, grid = g2)
mofa <- do_mofa(A, "meth_rna_17pairs", nf = 6)