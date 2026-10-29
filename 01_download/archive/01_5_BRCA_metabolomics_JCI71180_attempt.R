# =====================================================================
# STEP 1.5 — BREAST METABOLOMICS DOWNLOAD
# Output: brca_metabolites_raw.csv
#
# Source: Terunuma et al. 2014, J Clin Invest 124(1):398-412
#         (doi: 10.1172/JCI71180)
# 67 tumors + 65 tumor-adjacent noncancerous tissues, 536 metabolites
# (352 known, 184 unknown), Metabolon platform.
# Supplemental Tables 1-2 also carry clinical fields (incl. disease
# stage, ER status, race) needed later for early-stage subsetting.
# =====================================================================

# ---- 0. Setup ----
for (pkg in c("readxl", "httr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(readxl)
library(httr)

out_dir <- "data/breast_metabolomics"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)   # <- adjust to fit your project structure

# =====================================================================
# Download JCI supplemental metabolomics tables
# =====================================================================
cat("=== Downloading JCI71180 supplemental metabolomics tables ===\n")

jci_urls   <- c("https://www.jci.org/articles/view/71180/sd/1",
                "https://www.jci.org/articles/view/71180/sd/2")
dest_files <- c("JCI71180_supplement_1.xlsx",
                "JCI71180_supplement_2.xlsx")

ua <- "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0 Safari/537.36"

download_ok <- logical(length(jci_urls))
for (i in seq_along(jci_urls)) {
  cat("Attempting:", jci_urls[i], "\n")
  res <- tryCatch(
    httr::GET(jci_urls[i], httr::user_agent(ua),
              httr::write_disk(dest_files[i], overwrite = TRUE), httr::timeout(60)),
    error = function(e) e
  )
  download_ok[i] <- !inherits(res, "error") &&
    file.exists(dest_files[i]) &&
    file.size(dest_files[i]) > 1000
  cat(if (download_ok[i]) " -> downloaded OK\n" else " -> failed\n")
}

# =====================================================================
# Load, inspect, and save as the pipeline-standard raw file
# =====================================================================
if (any(download_ok)) {
  ok_file <- dest_files[download_ok][1]
  ext <- tools::file_ext(ok_file)

  if (ext %in% c("xlsx", "xls")) {
    sheets <- readxl::excel_sheets(ok_file)
    cat("Sheets found:", paste(sheets, collapse = ", "), "\n")
    metab_data <- readxl::read_excel(ok_file, sheet = sheets[1])

    # ---- Save as the Step 1.5 pipeline deliverable ----
    write.csv(metab_data, "brca_metabolites_raw.csv", row.names = FALSE)
    cat("Saved brca_metabolites_raw.csv (", nrow(metab_data), "rows x",
        ncol(metab_data), "cols )\n")

    # Quick look for a stage column, just for visibility - does NOT
    # filter the raw file. Subsetting/normalization happens in Step 2.4.
    stage_col <- grep("stage", colnames(metab_data), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(stage_col)) {
      cat("\nStage column detected:", stage_col, "\n")
      print(table(metab_data[[stage_col]]))
    } else {
      cat("\nNo obvious 'stage' column in this sheet - check other sheets in",
          ok_file, "or open it manually; tell me the exact header and I'll",
          "wire up early-stage filtering for Step 2.4.\n")
    }
  } else {
    cat("Downloaded file is not .xls/.xlsx (likely a PDF):", ok_file,
        "- open it manually to locate Supplemental Tables 1-2, save the",
        "metabolite table as brca_metabolites_raw.csv, then continue.\n")
  }
} else {
  message("
=========================== MANUAL STEP NEEDED ===========================
Automated download of the JCI supplemental tables failed - publisher
servers sometimes block scripted requests even with a browser User-Agent.

Manual fallback:
1. Go to: https://www.jci.org/articles/view/71180
2. Scroll to 'Supplemental material' near the bottom of the article
3. Download Supplemental Tables 1 and 2 (352 known + 184 unknown
   metabolites, 67 tumors / 65 adjacent noncancerous tissues)
4. Save into this folder, then run:
     metab_data <- readxl::read_excel('<saved_file>.xlsx')
     write.csv(metab_data, 'brca_metabolites_raw.csv', row.names = FALSE)
============================================================================
")
}

cat("\nStep 1.5 done. Output saved in:", getwd(), "\n")