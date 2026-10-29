## 05_0b  PRAD paired Normal-vs-Early cohort (Step 5.0)
## Early = Gleason score <= 7 (gleason_score in prad_rna_meta.rds); Normal = code 11; tumour = code 01
## Output: cohort_meth_rna.csv (17 pairs)
setwd("C:/Users/USER/Documents/MyThesis_MultiOmics")
out <- "05_FinalResults/Step5.0_Integration/PRAD"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
fix  <- function(x) gsub(".", "-", x, fixed = TRUE)
code <- function(b) substr(fix(b), 14, 15); pat <- function(b) substr(fix(b), 1, 12)

rm_  <- as.data.frame(readRDS("02_ProcessedData/Prostate/prad_rna_meta.rds"))
rna  <- data.frame(barcode = rm_$barcode, patient = rm_$patient, code = code(rm_$barcode))
early <- unique(rm_$patient[!is.na(rm_$gleason_score) & rm_$gleason_score <= 7])
cat("Patients with Gleason <=7:", length(early), "\n")

mh   <- fix(names(fread("02_ProcessedData/Prostate/prad_meth_normalized.csv", nrows = 0)))
mh   <- mh[grepl("^TCGA", mh)]
meth <- data.frame(barcode = mh, patient = pat(mh), code = code(mh))

both <- function(d) intersect(unique(d$patient[d$code == "11"]),
                              unique(d$patient[d$code == "01" & d$patient %in% early]))
P <- list(rna = both(rna), meth = both(meth))
print(sapply(P, length))
m2 <- intersect(P$meth, P$rna)
cat("PRAD Normal + Early (Gleason<=7) paired, meth+RNA:", length(m2), "\n")

pick <- function(d, p, g) sort(d$barcode[d$patient == p & d$code == ifelse(g == "Normal", "11", "01")])[1]
mk <- function(ps) do.call(rbind, lapply(ps, function(p) do.call(rbind, lapply(c("Normal","Early"), function(g)
  data.frame(patient = p, group = g, rna = pick(rna,p,g), meth = pick(meth,p,g))))))
coh <- mk(m2)
write.csv(coh, file.path(out, "cohort_meth_rna.csv"), row.names = FALSE)
print(coh)

