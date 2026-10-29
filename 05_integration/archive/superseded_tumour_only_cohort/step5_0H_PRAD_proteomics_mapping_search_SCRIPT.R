# ============================================================
# STEP 5.0H — FIND PRAD PROTEOMICS SAMPLE-ID MAPPING
# ============================================================

cat("\n============================================================\n")
cat("STEP 5.0H — SEARCHING PRAD PROTEOMICS FOR SAMPLE MAPPING\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 1. List ALL files in PRAD raw-data folder
# ------------------------------------------------------------

prad_files <- list.files(
  "01_RawData/Prostate",
  recursive = TRUE,
  full.names = TRUE
)

cat("\nTotal files found:", length(prad_files), "\n\n")

print(prad_files)


# ------------------------------------------------------------
# 2. Show files likely related to proteomics
# ------------------------------------------------------------

prot_files <- prad_files[
  grepl(
    "prot|protein|proteom|individual|sample|barcode|clinical|meta",
    basename(prad_files),
    ignore.case = TRUE
  )
]

cat("\n============================================================\n")
cat("POTENTIAL PROTEOMICS / SAMPLE-MAPPING FILES\n")
cat("============================================================\n")

print(prot_files)


# ------------------------------------------------------------
# 3. Inspect CSV files
# ------------------------------------------------------------

csv_files <- prad_files[
  grepl("\\.csv$", prad_files, ignore.case = TRUE)
]

cat("\n============================================================\n")
cat("CSV FILES\n")
cat("============================================================\n")

print(csv_files)


# ------------------------------------------------------------
# 4. Inspect dimensions and column names of every CSV
# ------------------------------------------------------------

if (length(csv_files) > 0) {

  for (f in csv_files) {

    cat("\n------------------------------------------------------------\n")
    cat("FILE:", f, "\n")
    cat("------------------------------------------------------------\n")

    x <- tryCatch(
      read.csv(
        f,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      error = function(e) NULL
    )

    if (!is.null(x)) {

      cat("Dimensions:", nrow(x), "x", ncol(x), "\n")

      cat("Columns:\n")
      print(colnames(x))

      cat("\nFirst 5 rows:\n")
      print(head(x, 5))

    } else {

      cat("Could not read as standard CSV.\n")

    }
  }
}


# ------------------------------------------------------------
# 5. Search text-readable files for "individual"
# ------------------------------------------------------------

text_files <- prad_files[
  grepl(
    "\\.(csv|txt|tsv|xlsx|xls)$",
    prad_files,
    ignore.case = TRUE
  )
]

cat("\n============================================================\n")
cat("SEARCHING FILE CONTENT FOR 'individual'\n")
cat("============================================================\n")

if (length(text_files) > 0) {

  for (f in text_files) {

    cat("\nChecking:", f, "\n")

    txt <- tryCatch(
      readLines(
        f,
        warn = FALSE,
        encoding = "UTF-8"
      ),
      error = function(e) character(0)
    )

    hits <- grep(
      "individual",
      txt,
      ignore.case = TRUE,
      value = TRUE
    )

    if (length(hits) > 0) {

      cat("MATCHES FOUND:\n")
      print(head(hits, 20))

    }

  }

}


# ------------------------------------------------------------
# 6. Search for TCGA identifiers
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SEARCHING FILE CONTENT FOR 'TCGA-'\n")
cat("============================================================\n")

if (length(text_files) > 0) {

  for (f in text_files) {

    txt <- tryCatch(
      readLines(
        f,
        warn = FALSE,
        encoding = "UTF-8"
      ),
      error = function(e) character(0)
    )

    hits <- grep(
      "TCGA-",
      txt,
      ignore.case = FALSE,
      value = TRUE
    )

    if (length(hits) > 0) {

      cat("\nFILE:", f, "\n")
      cat("Number of matching lines:", length(hits), "\n")
      print(head(hits, 10))

    }

  }

}


cat("\n============================================================\n")
cat("STEP 5.0H COMPLETE\n")
cat("============================================================\n")