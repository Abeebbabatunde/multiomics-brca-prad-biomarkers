## 05_2  Shuffled-label sanity check for DIABLO (Step 5.0)
## Run after 05_1a / 05_1b in the same R session (needs A, A2, A3, g2, g3, ntop, top_var).
## Normal/Early labels are swapped within a random half of the patients (seed 1); keepX from the real run.
## Reported errors (comp1 / comp2): BRCA 51 pairs 0.438 / 0.326; BRCA 16 pairs 0.303 / 0.178; PRAD 17 pairs 0.259 / 0.176.
suppressPackageStartupMessages(library(mixOmics))
sanity_check <- function(A, tag, folds, grid, kx_path) {
  set.seed(1)
  X <- mapply(top_var, A$X, ntop[names(A$X)], SIMPLIFY = FALSE)
  X <- lapply(X, function(m) withinVariation(t(m), design = data.frame(sample = A$patient)))
  design <- matrix(0.1, length(X), length(X), dimnames = list(names(X), names(X))); diag(design) <- 0
  kx <- readRDS(kx_path)$tune$choice.keepX
  Yp <- A$group; pts <- unique(A$patient)
  for (p in sample(pts, floor(length(pts) / 2))) { i <- which(A$patient == p); Yp[i] <- rev(Yp[i]) }
  m  <- block.splsda(X, Yp, ncomp = 2, keepX = kx, design = design)
  pf <- perf(m, validation = "Mfold", folds = folds, nrepeat = 10, dist = "centroids.dist")
  cat("\nSHUFFLED-LABEL sanity check -", tag, "\n")
  print(pf$WeightedVote.error.rate$centroids.dist)
  invisible(pf)
}
## PRAD (A from 05_1b)
sc_prad <- sanity_check(A, "meth_rna_17pairs", folds = 4, grid = g2,
                        kx_path = "05_FinalResults/Step5.0_Integration/PRAD/meth_rna_17pairs/DIABLO_objects.rds")
## BRCA (A2, A3 from 05_1a). Calls as run (sc2, sc3); the pasted function version for BRCA printed an extra
## line "Real-label error was 0..." which is cosmetic.
## sc2 <- sanity_check(A2, "meth_rna_51pairs", folds = 5, grid = g2, kx_path = file.path(base, "meth_rna_51pairs/DIABLO_objects.rds"))
## sc3 <- sanity_check(A3, "3omics_16pairs",  folds = 4, grid = g3, kx_path = file.path(base, "3omics_16pairs/DIABLO_objects.rds"))
