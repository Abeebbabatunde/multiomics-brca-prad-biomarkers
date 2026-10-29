## 05_0a  BRCA paired Normal-vs-Early cohorts (Step 5.0)
## Early = AJCC pathologic stage matching ^Stage (I|II)[AB]?$ ; Normal = sample code 11 ; tumour = code 01
## Outputs: cohort_3omics.csv (16 pairs) and cohort_meth_rna.csv (51 pairs)
setwd("C:/Users/USER/Documents/MyThesis_MultiOmics")
suppressPackageStartupMessages({ library(data.table) })
cat("\n## Cohort\n")
out <- "05_FinalResults/Step5.0_Integration/BRCA"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
fix  <- function(x) gsub(".", "-", x, fixed = TRUE)
code <- function(b) substr(fix(b), 14, 15); pat <- function(b) substr(fix(b), 1, 12)

rm_  <- as.data.frame(readRDS("02_ProcessedData/Breast/brca_rna_meta.rds"))
rna  <- data.frame(barcode = rm_$barcode, patient = rm_$patient, code = code(rm_$barcode))
early <- unique(rm_$patient[grepl("^Stage (I|II)[AB]?$", rm_$ajcc_pathologic_stage)])

mh   <- fix(names(fread("02_ProcessedData/Breast/brca_meth_normalized.csv", nrows = 0)))
mh   <- mh[grepl("^TCGA", mh)]
meth <- data.frame(barcode = mh, patient = pat(mh), code = code(mh))
mt   <- unique(meth$patient[meth$code == "01"])
cat("Stage of methylation tumor patients:\n")
print(table(rm_$ajcc_pathologic_stage[match(mt, rm_$patient)], useNA = "ifany"))

pl   <- fread("02_Normalized/Breast/brca_protein_normalized.csv")
pb   <- unique(fix(pl$sample_barcode))
prot <- data.frame(barcode = pb, patient = pat(pb), code = code(pb))

both <- function(d) intersect(unique(d$patient[d$code == "11"]),
                              unique(d$patient[d$code == "01" & d$patient %in% early]))
P <- list(rna = both(rna), meth = both(meth), prot = both(prot))
print(sapply(P, length))
all3 <- Reduce(intersect, P); m2 <- intersect(P$meth, P$rna)
cat("3 omics (paired Normal + Early):", length(all3), "| meth+RNA only:", length(m2), "\n")

pick <- function(d, p, g) sort(d$barcode[d$patient == p & d$code == ifelse(g == "Normal", "11", "01")])[1]
mk <- function(ps, withprot) do.call(rbind, lapply(ps, function(p) do.call(rbind, lapply(c("Normal","Early"), function(g)
  data.frame(patient = p, group = g, rna = pick(rna,p,g), meth = pick(meth,p,g),
             prot = if (withprot) pick(prot,p,g) else NA)))))
write.csv(mk(all3, TRUE),  file.path(out, "cohort_3omics.csv"),   row.names = FALSE)
write.csv(mk(m2,  FALSE),  file.path(out, "cohort_meth_rna.csv"), row.names = FALSE)
