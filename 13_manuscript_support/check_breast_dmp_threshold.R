# Recount the breast differentially methylated positions at both delta-beta thresholds.
# No model is re-run: the counts come from the saved limma table of Step 3.2.
# In that table logFC is the difference in mean beta (tumour minus normal).

project_dir <- "C:/Users/USER/Documents/MyThesis_MultiOmics"   # edit if different

f <- list.files(project_dir, pattern = "breast_methylation_diffmeth_results.csv",
                recursive = TRUE, full.names = TRUE)[1]
stopifnot(!is.na(f))
res <- read.csv(f)
stopifnot(all(c("logFC", "adj.P.Val") %in% names(res)))

count_dmp <- function(cut) {
  keep <- !is.na(res$adj.P.Val) & res$adj.P.Val < 0.05 & abs(res$logFC) >= cut
  c(threshold = cut, total = sum(keep),
    hyper = sum(keep & res$logFC > 0), hypo = sum(keep & res$logFC < 0))
}

out <- rbind(count_dmp(0.1), count_dmp(0.2))
cat("Probes tested:", nrow(res), "\n")
cat("Adjusted P < 0.05 only:", sum(res$adj.P.Val < 0.05, na.rm = TRUE), "\n\n")
print(out)
write.csv(out, "breast_DMP_threshold_sensitivity.csv", row.names = FALSE)
