## 05_4a  BRCA gene-level overlap: methylation (DMP genes), transcriptome (DEGs), proteome (significant DEPs)
## Results reported: 4,763 / 7,924 / 109 genes; pairwise 1,509 (meth-RNA), 19 (meth-prot), 29 (RNA-prot); 8 genes in all three layers.
## SETUP LINES: taken from the script as run in the console session (first attempt stopped at the missing
##   UpSetR / ggVennDiagram / annotation-package errors; the packages were then installed and the minfi/annotation
##   packages attached with the library() line below, after which the rest ran and produced the reported counts).
suppressPackageStartupMessages({ library(data.table); library(UpSetR); library(ggVennDiagram); library(ggplot2); library(org.Hs.eg.db) })
for (p in c("minfi","IlluminaHumanMethylation450kanno.ilmn12.hg19"))
  if (!requireNamespace(p, quietly = TRUE)) stop("Install first: BiocManager::install('", p, "')")
setwd("C:/Users/USER/Documents/MyThesis_MultiOmics")
o <- "05_FinalResults/Step5.0_Integration/BRCA"; dir.create(o, recursive = TRUE, showWarnings = FALSE)
sets <- list()

suppressPackageStartupMessages({
  library(minfi)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
})

dmp <- fread("03_Results/Epigenomics/brca_DMPs.csv")
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
g   <- ann[match(dmp$Probe, rownames(ann)), "UCSC_RefGene_Name"]
cat("DMP probes found in annotation:", sum(dmp$Probe %in% rownames(ann)), "of", nrow(dmp), "\n")
sets$Methylation <- unique(na.omit(trimws(unlist(strsplit(as.character(g), ";")))))
deg <- fread("03_Results/Transcriptomics/brca_DEGs.csv"); print(head(deg$V1, 3))
sets$Transcriptome <- unique(na.omit(mapIds(org.Hs.eg.db, sub("\\..*", "", deg$V1), "SYMBOL", "ENSEMBL", multiVals = "first")))

dep <- fread("03_Results/Proteomics/brca_DEPs.csv"); sg <- dep$significant
dep <- dep[if (is.logical(sg)) which(sg) else which(tolower(as.character(sg)) %in% c("true","yes","significant","sig","1")), ]
stopifnot(nrow(dep) > 0); cat("Significant DEPs:", nrow(dep), "\n")
cl <- toupper(sub("_p[STY].*$", "", dep$peptide_target))
alias <- c("1433BETA"="YWHAB","1433EPSILON"="YWHAE","1433ZETA"="YWHAZ","4EBP1"="EIF4EBP1","ACC1"="ACACA",
  "ERALPHA"="ESR1","HER2"="ERBB2","PR"="PGR","P53"="TP53","P27"="CDKN1B","ECADHERIN"="CDH1",
  "PCADHERIN"="CDH3","BETACATENIN"="CTNNB1","CYCLINE1"="CCNE1","CYCLINB1"="CCNB1",
  "B7-H3"="CD276","CAVEOLIN1"="CAV1","CD134"="TNFRSF4","MCL-1"="MCL1","MITOFUSIN-1"="MFN1",
  "MITOFUSIN-2"="MFN2","HEXOKINASE-I"="HK1","HEXOKINASE-II"="HK2","NOTCH1-CLEAVED"="NOTCH1",
  "PERK"="EIF2AK3","ATP5H"="ATP5PD","GSK-3B"="GSK3B","GSK3ALPHABETA"="GSK3A;GSK3B","PORIN"="VDAC1",
  "HMHA1"="ARHGAP45","CMYC"="MYC","COX-IV"="COX4I1","DNA_POLG"="POLG","DRP1"="DNM1L","HSP27"="HSPB1",
  "HSP70"="HSPA1A","RIP"="RIPK1","EIF4G"="EIF4G1","IRF-3"="IRF3","RSK1"="RPS6KA1","EMA"="MUC1",
  "SHP2"="PTPN11","RPA32"="RPA2","STING"="STING1","ZAP-70"="ZAP70","P70S6K"="RPS6KB1",
  "COUP-TFII"="NR2F2","MCT4"="SLC16A3","MNK1"="MKNK1","IL-6"="IL6","FAK"="PTK2","COLLAGENVI"="COL6A1",
  "MERIT40"="BABAM1","FGF-BASIC"="FGF2","PLC-GAMMA1"="PLCG1","CD171"="L1CAM","CALNEXIN"="CANX",
  "53BP1"="TP53BP1","CYCLOPHILIN-F"="PPIF","S6"="RPS6","SMAC"="DIABLO","CLAUDIN7"="CLDN7",
  "GLUTAMATE-D1-2"="GLUD1","CREB"="CREB1","MTCO1"="MT-CO1")
cl  <- ifelse(cl %in% names(alias), alias[cl], cl)
gs  <- unlist(strsplit(cl, ";"))
ok  <- gs %in% keys(org.Hs.eg.db, "SYMBOL")
cat("Mapped:", sum(cl %in% names(alias) | cl %in% keys(org.Hs.eg.db, "SYMBOL")), "of", length(cl), "\n")
cat("Still unmatched:", unique(gs[!ok]), "\n")
sets$Proteome <- unique(gs[ok])

print(sapply(sets, length))
cmb <- combn(names(sets), 2, simplify = FALSE)
print(sapply(cmb, function(p) setNames(length(intersect(sets[[p[1]]], sets[[p[2]]])), paste(p, collapse = " & "))))
cat("In all layers:", length(Reduce(intersect, sets)), "\n")

ggsave(file.path(o, "venn.png"), ggVennDiagram(sets, label = "count") +
       scale_fill_gradient(low = "white", high = "steelblue"), width = 6, height = 5, dpi = 300)
png(file.path(o, "upset.png"), 2000, 1400, res = 300)
print(upset(fromList(sets), order.by = "freq", nsets = length(sets), text.scale = 1.3)); dev.off()
write.csv(data.frame(gene = Reduce(intersect, sets)), file.path(o, "genes_in_all_layers.csv"), row.names = FALSE)
