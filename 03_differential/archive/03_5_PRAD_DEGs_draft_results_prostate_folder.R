## Earlier draft of Step 3.5 (writes to results/prostate/). Its outputs are not in the project listing, so it did not produce the files used later.
library(DESeq2); library(dplyr); library(tibble)
dds_prad <- readRDS("02_ProcessedData/Prostate/dds_prad.rds")
class(dds_prad); design(dds_prad); colData(dds_prad); resultsNames(dds_prad)   # printed: design ~group, 219 samples, resultsNames = character(0)
if (length(resultsNames(dds_prad)) == 0) {
  dds_prad <- DESeq(dds_prad)
  saveRDS(dds_prad, "02_ProcessedData/Prostate/dds_prad.rds")
}
res <- results(dds_prad, contrast = c("group", "EarlyProstateCancer", "Normal"), alpha = 0.05)
res <- lfcShrink(dds_prad, contrast = c("group", "EarlyProstateCancer", "Normal"), res = res, type = "ashr")
res_df <- as.data.frame(res) %>% rownames_to_column("ensembl_gene_id") %>% arrange(padj)
prad_DEGs <- res_df %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1)
dir.create("results/prostate", recursive = TRUE, showWarnings = FALSE)
write.csv(res_df, "results/prostate/PRAD_DESeq2_all_results.csv", row.names = FALSE)
write.csv(prad_DEGs, "results/prostate/PRAD_significant_DEGs.csv", row.names = FALSE)
