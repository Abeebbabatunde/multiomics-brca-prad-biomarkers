## 05_4b  PRAD gene-level overlap, methylation (DMP genes) and transcriptome (DEGs)
## Uses prad_DMPs.csv (63,993 DMPs). genes_in_all_layers.csv holds the TWO-layer overlap for PRAD.
setwd("C:/Users/USER/Documents/MyThesis_MultiOmics")
suppressPackageStartupMessages({ library(data.table); library(UpSetR); library(ggVennDiagram); library(ggplot2)
  library(org.Hs.eg.db); library(minfi); library(IlluminaHumanMethylation450kanno.ilmn12.hg19) })
o <- "05_FinalResults/Step5.0_Integration/PRAD"; dir.create(o, recursive = TRUE, showWarnings = FALSE)
sets <- list()

dmp <- fread("03_Results/Epigenomics/prad_DMPs.csv")
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
g   <- ann[match(dmp$probe_id, rownames(ann)), "UCSC_RefGene_Name"]
cat("DMP probes found in annotation:", sum(dmp$probe_id %in% rownames(ann)), "of", nrow(dmp), "\n")
sets$Methylation <- unique(na.omit(trimws(unlist(strsplit(as.character(g), ";")))))

deg <- fread("03_Results/Transcriptomics/prad_DEGs.csv")
sets$Transcriptome <- unique(na.omit(mapIds(org.Hs.eg.db, sub("\\..*", "", deg$ensembl_gene_id),
                                            "SYMBOL", "ENSEMBL", multiVals = "first")))

print(sapply(sets, length))
cat("Overlap (both layers):", length(intersect(sets$Methylation, sets$Transcriptome)), "\n")

ggsave(file.path(o, "venn.png"), ggVennDiagram(sets, label = "count") +
         scale_fill_gradient(low = "white", high = "steelblue"), width = 6, height = 5, dpi = 300)
png(file.path(o, "upset.png"), 2000, 1400, res = 300)
print(upset(fromList(sets), order.by = "freq", nsets = length(sets), text.scale = 1.3)); dev.off()
write.csv(data.frame(gene = intersect(sets$Methylation, sets$Transcriptome)),
          file.path(o, "genes_in_all_layers.csv"), row.names = FALSE)
