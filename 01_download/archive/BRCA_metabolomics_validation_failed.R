# ================================================
# STEP 1.5 COMPLETE VALIDATION
# ================================================

cat("========================================\n")
cat("   STEP 1.5 VALIDATION\n")
cat("========================================\n\n")

# -----------------------------------------------
# CHECK 1: File exists and size
# -----------------------------------------------
cat("CHECK 1: File integrity...\n\n")

f1 <- "01_RawData/Breast/brca_metabolites_raw.csv"
f2 <- "01_RawData/Breast/brca_metabolites_meta.csv"

cat(sprintf("  Data file  : %s  %.2f MB\n",
            ifelse(file.exists(f1),"✅","❌"),
            ifelse(file.exists(f1),
              file.size(f1)/1024/1024, 0)))
cat(sprintf("  Meta file  : %s  %.2f MB\n\n",
            ifelse(file.exists(f2),"✅","❌"),
            ifelse(file.exists(f2),
              file.size(f2)/1024/1024, 0)))

if(!file.exists(f1)){
  cat("❌ Data file missing — re-run Step 1.5\n")
  stop("File not found")
}

# -----------------------------------------------
# CHECK 2: Load and inspect
# -----------------------------------------------
cat("CHECK 2: Loading data...\n\n")

brca_metabolites_early <- read.csv(
  f1, row.names=1, check.names=FALSE)
brca_met_meta <- read.csv(f2,
  check.names=FALSE)

cat(sprintf("  Metabolites (rows): %d\n",
            nrow(brca_metabolites_early)))
cat(sprintf("  Samples (cols)    : %d\n\n",
            ncol(brca_metabolites_early)))

# -----------------------------------------------
# CHECK 3: Missing values
# -----------------------------------------------
cat("CHECK 3: Missing values...\n\n")

missing_pct <- round(
  100 * sum(is.na(brca_metabolites_early)) /
  (nrow(brca_metabolites_early) *
   ncol(brca_metabolites_early)), 1)

cat(sprintf("  Missing: %.1f%%  %s\n\n",
            missing_pct,
            ifelse(missing_pct < 50,
                   "✅ Acceptable",
                   "❌ Too high")))

# -----------------------------------------------
# CHECK 4: Show actual values
# -----------------------------------------------
cat("CHECK 4: Actual data values...\n\n")

cat("First 5 metabolites, first 4 samples:\n")
print(brca_metabolites_early[
  1:min(5, nrow(brca_metabolites_early)),
  1:min(4, ncol(brca_metabolites_early))])

cat("\nValue range:\n")
val_range <- range(brca_metabolites_early,
                    na.rm=TRUE)
cat(sprintf("  Min: %.4f\n", val_range[1]))
cat(sprintf("  Max: %.4f\n\n", val_range[2]))

# -----------------------------------------------
# CHECK 5: Group distribution
# -----------------------------------------------
cat("CHECK 5: Group distribution...\n\n")

cat(sprintf("  Cancer (MB_) : %d\n",
  sum(brca_met_meta$Group ==
      "EarlyBreastCancer")))
cat(sprintf("  Normal (HCB_): %d\n\n",
  sum(brca_met_meta$Group == "Normal")))

# -----------------------------------------------
# CHECK 6: Sample names format
# -----------------------------------------------
cat("CHECK 6: Sample names...\n\n")

cat("First 8 sample names:\n")
print(head(colnames(brca_metabolites_early), 8))

# -----------------------------------------------
# VERDICT
# -----------------------------------------------
cat("\n========================================\n")
cat("   VALIDATION VERDICT\n")
cat("========================================\n\n")

checks <- list(
  list("File exists",
       file.exists(f1)),
  list("Has metabolites (>0)",
       nrow(brca_metabolites_early) > 0),
  list("Has samples (>10)",
       ncol(brca_metabolites_early) > 10),
  list("Has cancer samples",
       sum(brca_met_meta$Group ==
           "EarlyBreastCancer") > 0),
  list("Has normal samples",
       sum(brca_met_meta$Group ==
           "Normal") > 0),
  list("Missing values acceptable (<50%)",
       missing_pct < 50),
  list("Has real values (not all NA)",
       sum(!is.na(brca_metabolites_early)) > 0)
)

all_pass <- TRUE
for(chk in checks){
  pass <- chk[[2]]
  cat(sprintf("  %-40s %s\n",
              chk[[1]],
              ifelse(pass,
                     "✅ PASS",
                     "❌ FAIL")))
  if(!pass) all_pass <- FALSE
}

cat("\n")
if(all_pass){
  cat("  ✅ ALL CHECKS PASSED\n")
  cat("  Data is good for analysis\n")
  cat("  Continue to Phase 2\n")
} else {
  cat("  ❌ SOME CHECKS FAILED\n")
  cat("  Tell me which checks failed\n")
}
cat("========================================\n")