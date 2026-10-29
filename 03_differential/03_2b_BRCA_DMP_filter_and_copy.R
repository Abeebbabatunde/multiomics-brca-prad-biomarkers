# STEP 3.2b — BRCA DMP filter on the Step 3.2 limma table (AS RUN, console: nrow(sig_results) = 18565).
# Requires diff_results and output_dir from 03_2_BRCA_methylation_limma_AS_RUN.R in the same session.
delta_beta_cutoff <- 0.2   # (comment in original: "try 0.1 if you want a larger candidate list")
sig_results <- diff_results[diff_results$adj.P.Val < 0.05 & abs(diff_results$logFC) > delta_beta_cutoff, ]
nrow(sig_results)          # 18565
write.csv(sig_results, file.path(output_dir, "breast_methylation_DMPs_filtered.csv"), row.names = FALSE)

# Copy as pasted (the copy command that was pasted):
dir.create("03_DifferentialAnalysis/Breast", recursive = TRUE, showWarnings = FALSE)
file.copy("results/step2_2_methylation_diff/breast_methylation_DMPs_filtered.csv",
          "03_DifferentialAnalysis/Breast/brca_DMPs.csv")
# NOTE: Step 5 prep scripts read "03_Results/Epigenomics/brca_DMPs.csv". The command that put the file there was NOT seen
# (the Phase-3 template Step 3.2 writes to that path, but that template is not the script that produced these results).
# Confirm that 03_Results/Epigenomics/brca_DMPs.csv is a copy of this filtered file (18,565 rows).
