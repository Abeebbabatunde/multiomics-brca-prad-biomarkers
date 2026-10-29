## 05_3  DIABLO signatures annotated with the per-omics differential results (BRCA)
## Output: signature_with_DE.csv for meth_rna_51pairs and 3omics_16pairs
## Reported: 51 pairs meth 9/20 and RNA 10/20 significant; 16 pairs meth 2/20, prot 7/10, RNA 15/20.
## NOTE: de_effect/de_padj are NA for features absent from the significant-only DMP table.
suppressPackageStartupMessages({ library(data.table); library(org.Hs.eg.db) })
setwd("C:/Users/USER/Documents/MyThesis_MultiOmics")
b  <- "05_FinalResults/Step5.0_Integration/BRCA"
s2 <- fread(file.path(b, "meth_rna_51pairs/DIABLO_selected_features.csv"))
s3 <- fread(file.path(b, "3omics_16pairs/DIABLO_selected_features.csv"))

deg <- fread("03_Results/Transcriptomics/brca_DEGs.csv")
deg$symbol <- mapIds(org.Hs.eg.db, sub("\\..*", "", deg$V1), "SYMBOL", "ENSEMBL", multiVals = "first")
dmp <- fread("03_Results/Epigenomics/brca_DMPs.csv")
dep <- fread("03_Results/Proteomics/brca_DEPs.csv")

annot <- function(s) {
  s$de_effect <- NA_real_; s$de_padj <- NA_real_
  i <- s$block == "rna";  k <- match(s$feature[i], deg$symbol);         s$de_effect[i] <- deg$log2FoldChange[k]; s$de_padj[i] <- deg$padj[k]
  i <- s$block == "meth"; k <- match(s$feature[i], dmp$Probe);          s$de_effect[i] <- dmp$logFC[k];          s$de_padj[i] <- dmp$adj.P.Val[k]
  i <- s$block == "prot"; k <- match(s$feature[i], dep$peptide_target); s$de_effect[i] <- dep$mean_diff[k];      s$de_padj[i] <- dep$p_adj[k]
  s$sig_in_DE <- !is.na(s$de_padj) & s$de_padj < 0.05
  s[order(block, comp, -abs(loading))]
}
a2 <- annot(s2); a3 <- annot(s3)
print(as.data.frame(a2), row.names = FALSE)
cat("\n51-pair: features also significant in your DE steps:\n"); print(table(a2$block, a2$sig_in_DE))
cat("\n16-pair:\n"); print(table(a3$block, a3$sig_in_DE))
fwrite(a2, file.path(b, "meth_rna_51pairs/signature_with_DE.csv"))
fwrite(a3, file.path(b, "3omics_16pairs/signature_with_DE.csv"))
