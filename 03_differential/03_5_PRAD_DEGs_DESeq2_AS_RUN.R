## ================================================
## STEP 3.5 — Prostate DEGs (TCGA-PRAD)
## Output: 03_Results/Transcriptomics/prad_DEGs.csv
## ================================================

setwd("MyThesis_MultiOmics")
source("00_Scripts/00_STARTUP.R")   # reloads dds_prad, prad_rna, etc.

cat("\n========================================\n")
cat("   STEP 3.5 — PROSTATE DEGs\n")
cat("========================================\n\n")

## ---- 1. Load dds_prad (already has size factors estimated) ----
dds_prad <- readRDS("02_ProcessedData/Prostate/dds_prad.rds")

design(dds_prad)          # ~group
table(colData(dds_prad)$group)   # EarlyProstateCancer vs Normal

## ---- 2. Run DESeq() if not already fitted ----
if (length(resultsNames(dds_prad)) == 0) {
  cat("Fitting DESeq2 model (dispersions + GLM)...\n")
  dds_prad <- DESeq(dds_prad)
  saveRDS(dds_prad, "02_ProcessedData/Prostate/dds_prad.rds")
  cat("  ✅ dds_prad.rds updated with fitted model\n")
} else {
  cat("  ✅ DESeq() already fitted —", resultsNames(dds_prad), "\n")
}

## ---- 3. Extract results: EarlyProstateCancer vs Normal ----
res <- results(dds_prad,
                contrast = c("group", "EarlyProstateCancer", "Normal"),
                alpha = 0.05)

res <- lfcShrink(dds_prad,
                  contrast = c("group", "EarlyProstateCancer", "Normal"),
                  res = res, type = "ashr")

prad_res_full <- as.data.frame(res) %>%
  rownames_to_column("ensembl_gene_id") %>%
  arrange(padj)

## ---- 4. Define significant DEGs ----
padj_cutoff <- 0.05
lfc_cutoff  <- 1

prad_DEGs <- prad_res_full %>%
  filter(!is.na(padj), padj < padj_cutoff, abs(log2FoldChange) >= lfc_cutoff)

cat("\nTotal significant Prostate DEGs:", nrow(prad_DEGs), "\n")
cat("  Upregulated (EarlyProstateCancer):", sum(prad_DEGs$log2FoldChange > 0), "\n")
cat("  Downregulated (EarlyProstateCancer):", sum(prad_DEGs$log2FoldChange < 0), "\n")

## ---- 5. Save outputs (matches your step tracker path) ----
dir.create("03_Results/Transcriptomics", recursive = TRUE, showWarnings = FALSE)

write.csv(prad_res_full,
          "03_Results/Transcriptomics/prad_DEGs_full_results.csv",
          row.names = FALSE)

write.csv(prad_DEGs,
          "03_Results/Transcriptomics/prad_DEGs.csv",   # <- matches 00_STARTUP.R tracker check
          row.names = FALSE)

cat("\n  ✅ Saved: 03_Results/Transcriptomics/prad_DEGs.csv\n")

## ---- 6. Optional QC volcano plot ----
## (needs BiocManager::install("EnhancedVolcano") + library(EnhancedVolcano); the first attempt failed with "could not find function")
dir.create("04_Figures", recursive = TRUE, showWarnings = FALSE)
png("04_Figures/prad_DEG_volcano.png", width = 1000, height = 800, res = 120)
EnhancedVolcano(prad_res_full,
  lab = prad_res_full$ensembl_gene_id,
  x = 'log2FoldChange',
  y = 'padj',
  pCutoff = padj_cutoff,
  FCcutoff = lfc_cutoff,
  title = 'TCGA-PRAD: EarlyProstateCancer vs Normal',
  subtitle = NULL)
dev.off()
cat("  ✅ Saved: 04_Figures/prad_DEG_volcano.png\n")

cat("\n========================================\n")
cat("   STEP 3.5 COMPLETE\n")
cat("   ▶️  NEXT: Step 3.6 — Prostate DMPs\n")
cat("========================================\n")
