# STEP 10.3b: RNA log2 fold change, CpG annotation, pathway overlap (run after 10_1_final_panel_tables.R)
project_root <- "C:/Users/USER/Documents/MyThesis_MultiOmics"
P <- function(...) file.path(project_root, ...)
out_dir <- P("10_FinalResults", "Step10_FinalPanel")
step7 <- readRDS(P("07_FinalResults","Step7_FeatureSelection_ML","Step7_feature_selection_results.rds"))
panel <- list(BRCA = step7$BRCA$consensus, PRAD = step7$PRAD$consensus)

# (1) RNA paired log2 fold change
l2fc <- function(cn) {
  g <- panel[[cn]][!grepl("^cg", panel[[cn]]) & panel[[cn]] != "VAV1"]
  d <- as.data.frame(readRDS(P("07_FinalResults","Step7_FeatureSelection_ML", paste0(cn, "_targeted_dataset.rds"))))
  data.frame(cancer = cn, gene = g, median_log2FC = sapply(g, function(x) {
    m <- merge(d[d$group == "Normal", c("patient", x)], d[d$group == "Early", c("patient", x)], by = "patient")
    round(median(log2((m[[3]] + 1) / (m[[2]] + 1)), na.rm = TRUE), 2) }), row.names = NULL)
}
fc <- rbind(l2fc("BRCA"), l2fc("PRAD")); print(fc, row.names = FALSE)
write.csv(fc, file.path(out_dir, "Step10_RNA_log2FC.csv"), row.names = FALSE)

# (2) CpG -> gene annotation (needs the 450K annotation package)
pk <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
if (requireNamespace(pk, quietly = TRUE)) {
  e <- new.env(); data("Other", package = pk, envir = e)
  cp <- unique(unlist(lapply(panel, function(x) x[grepl("^cg", x)]), use.names = FALSE))
  an <- as.data.frame(e$Other[cp, c("UCSC_RefGene_Name", "UCSC_RefGene_Group")])
  an$probe <- rownames(an); print(an, row.names = FALSE)
  write.csv(an, file.path(out_dir, "Step10_CpG_annotation.csv"), row.names = FALSE)
  gs <- function(cn) {
    p <- intersect(panel[[cn]], an$probe)
    g <- unique(unlist(strsplit(as.character(an$UCSC_RefGene_Name[match(p, an$probe)]), ";")))
    unique(c(g[nzchar(g)], panel[[cn]][!grepl("^cg", panel[[cn]])]))
  }
  cat("\nShared genes after CpG annotation:", paste(intersect(gs("BRCA"), gs("PRAD")), collapse = ", "), "\n")
} else cat("\nCpG annotation skipped: install", pk, "via BiocManager\n")

# (3) pathway overlap from Step 6 enrichment (FDR < 0.05)
ff <- list.files(P("06_FinalResults", "Step6_Network"), pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
ff <- ff[grepl("GO_BP|KEGG|Hallmark|Reactome", ff) & !grepl("top20|overlap", ff, ignore.case = TRUE)]
cat("\nEnrichment files:\n"); print(sub(P(""), "", ff, fixed = TRUE))
read_enr <- function(f) {
  x <- tryCatch(read.csv(f, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(x) || !nrow(x)) return(NULL)
  tc <- intersect(c("Description", "Term", "pathway", "ID"), names(x))[1]
  pc <- intersect(c("p.adjust", "padj", "FDR", "qvalue", "adj.P.Val"), names(x))[1]
  if (is.na(tc) || is.na(pc)) { cat("  cannot read:", basename(f), "| cols:", paste(names(x), collapse = ", "), "\n"); return(NULL) }
  x <- x[!is.na(x[[pc]]) & x[[pc]] < 0.05, ]
  if (!nrow(x)) return(NULL)
  data.frame(cancer = ifelse(grepl("BRCA", f), "BRCA", ifelse(grepl("PRAD", f), "PRAD", NA)),
             db = sub(".*(GO_BP|KEGG|Hallmark|Reactome).*", "\\1", basename(f)),
             term = tolower(trimws(x[[tc]])), stringsAsFactors = FALSE)
}
enr <- do.call(rbind, lapply(ff, read_enr))
if (!is.null(enr)) {
  enr <- unique(enr[!is.na(enr$cancer), ])
  ov <- do.call(rbind, lapply(unique(enr$db), function(d) {
    b <- enr$term[enr$cancer == "BRCA" & enr$db == d]; p <- enr$term[enr$cancer == "PRAD" & enr$db == d]
    data.frame(db = d, n_BRCA = length(b), n_PRAD = length(p), n_shared = length(intersect(b, p)),
               shared = paste(head(intersect(b, p), 15), collapse = " | "))
  }))
  print(ov, row.names = FALSE)
  write.csv(ov, file.path(out_dir, "Step10_pathway_overlap.csv"), row.names = FALSE)
} else cat("No significant enrichment terms could be read.\n")
