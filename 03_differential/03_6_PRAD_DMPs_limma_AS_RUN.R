## ================================================
## STEP 3.6 — Prostate DMPs (TCGA-PRAD methylation)
## Output: 03_Results/Epigenomics/prad_DMPs.csv
## Method: chunked limma on M-values (same approach as Breast DMPs)
## ================================================

setwd("MyThesis_MultiOmics")
source("00_Scripts/00_STARTUP.R")

cat("\n========================================\n")
cat("   STEP 3.6 — PROSTATE DMPs\n")
cat("========================================\n\n")

## ---- 1. Load methylation object + metadata ----
prad_meth        <- readRDS("01_RawData/Prostate/prad_meth.rds")
prad_meth_meta   <- readRDS("02_ProcessedData/Prostate/prad_meth_meta.rds")

## ---- 2. Inspect structure BEFORE running (confirm these match) ----
class(prad_meth)
colnames(prad_meth_meta)   # <- confirm group column name here
table(prad_meth_meta$group)  # <- adjust column name below if different

## ---- 3. Extract beta-value matrix ----
if (is(prad_meth, "SummarizedExperiment") || is(prad_meth, "GenomicRatioSet")) {
  beta_matrix <- assay(prad_meth)
} else {
  beta_matrix <- as.matrix(prad_meth)
}

# Ensure metadata rows align with matrix columns
prad_meth_meta <- prad_meth_meta[match(colnames(beta_matrix), rownames(prad_meth_meta)), ]
stopifnot(all(colnames(beta_matrix) == rownames(prad_meth_meta)))

cat("Methylation matrix:", nrow(beta_matrix), "probes x",
    ncol(beta_matrix), "samples\n")

## ---- 4. Convert beta -> M-values (limma performs better on M-values) ----
# Clip betas away from 0/1 to avoid -Inf/Inf
beta_matrix[beta_matrix <= 0] <- 1e-6
beta_matrix[beta_matrix >= 1] <- 1 - 1e-6
m_matrix <- log2(beta_matrix / (1 - beta_matrix))

## ---- 5. Design matrix ----
group <- factor(prad_meth_meta$group, levels = c("Normal", "EarlyProstateCancer"))
design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "EarlyProstateCancer_vs_Normal")

## ---- 6. Chunked limma (memory-safe for large probe counts) ----
chunk_size <- 50000
n_probes   <- nrow(m_matrix)
chunks     <- split(seq_len(n_probes), ceiling(seq_len(n_probes) / chunk_size))

cat("Running limma in", length(chunks), "chunks of ~", chunk_size, "probes...\n")

dmp_results_list <- vector("list", length(chunks))

for (i in seq_along(chunks)) {
  idx <- chunks[[i]]
  fit <- lmFit(m_matrix[idx, , drop = FALSE], design)
  fit <- eBayes(fit)
  tt  <- topTable(fit, coef = "EarlyProstateCancer_vs_Normal",
                   number = length(idx), sort.by = "none")
  tt$probe_id <- rownames(m_matrix)[idx]

  # also compute mean delta-beta (effect size in beta units, easier to interpret)
  beta_early <- rowMeans(beta_matrix[idx, group == "EarlyProstateCancer", drop = FALSE])
  beta_norm  <- rowMeans(beta_matrix[idx, group == "Normal", drop = FALSE])
  tt$delta_beta <- beta_early - beta_norm

  dmp_results_list[[i]] <- tt
  cat("  Chunk", i, "/", length(chunks), "done\n")
}

prad_dmp_full <- do.call(rbind, dmp_results_list) %>%
  select(probe_id, logFC, delta_beta, P.Value, adj.P.Val, everything()) %>%
  arrange(adj.P.Val)

## ---- 7. Define significant DMPs ----
padj_cutoff       <- 0.05
delta_beta_cutoff <- 0.1   # standard epigenetics effect-size threshold

prad_DMPs <- prad_dmp_full %>%
  filter(adj.P.Val < padj_cutoff, abs(delta_beta) >= delta_beta_cutoff)

cat("\nTotal significant Prostate DMPs:", nrow(prad_DMPs), "\n")
cat("  Hypermethylated (EarlyProstateCancer):", sum(prad_DMPs$delta_beta > 0), "\n")
cat("  Hypomethylated (EarlyProstateCancer):", sum(prad_DMPs$delta_beta < 0), "\n")

## ---- 8. Save outputs (matches step tracker path) ----
dir.create("03_Results/Epigenomics", recursive = TRUE, showWarnings = FALSE)

write.csv(prad_dmp_full,
          "03_Results/Epigenomics/prad_DMPs_full_results.csv",
          row.names = FALSE)

write.csv(prad_DMPs,
          "03_Results/Epigenomics/prad_DMPs.csv",   # <- matches 00_STARTUP.R tracker check
          row.names = FALSE)

cat("\n  ✅ Saved: 03_Results/Epigenomics/prad_DMPs.csv\n")

cat("\n========================================\n")
cat("   STEP 3.6 COMPLETE\n")
cat("   ▶️  NEXT: Step 3.7 — Prostate DEPs\n")
cat("========================================\n")