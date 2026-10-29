# ============================================================
# STEP 5.0I — FIND PRAD PROTEOMICS → TCGA MAPPING
# ============================================================

cat("\n============================================================\n")
cat(" SEARCHING FOR PRAD PROTEOMICS TCGA MAPPING\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Search filenames throughout project
# ------------------------------------------------------------

all_files <- list.files(
  ".",
  recursive = TRUE,
  full.names = TRUE
)

candidate_files <- all_files[
  grepl(
    "prad|prostate|proteom|clinical|metadata|sample|patient|mapping|annotation|E-PROT|TCGA",
    basename(all_files),
    ignore.case = TRUE
  )
]

cat("\nCandidate files:\n")
print(candidate_files)


# ------------------------------------------------------------
# 2. Search text-based metadata files for "individual"
#    and "TCGA"
# ------------------------------------------------------------

text_files <- candidate_files[
  grepl(
    "\\.(csv|tsv|txt|xlsx|xls)$",
    candidate_files,
    ignore.case = TRUE
  )
]

cat("\n============================================================\n")
cat(" SEARCHING METADATA FILES\n")
cat("============================================================\n")

for (f in text_files) {

  cat("\n------------------------------------------------------------\n")
  cat("FILE:", f, "\n")
  cat("------------------------------------------------------------\n")

  try({

    if (grepl("\\.csv$", f, ignore.case = TRUE)) {

      dat <- read.csv(
        f,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

    } else if (grepl("\\.tsv$", f, ignore.case = TRUE)) {

      dat <- read.delim(
        f,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

    } else {

      next
    }

    # Show dimensions
    cat(
      "Dimensions:",
      nrow(dat),
      "rows ×",
      ncol(dat),
      "columns\n"
    )

    # Show potentially useful columns
    useful_cols <- colnames(dat)[
      grepl(
        "individual|sample|patient|case|barcode|submitter|TCGA|subject|id",
        colnames(dat),
        ignore.case = TRUE
      )
    ]

    if (length(useful_cols) > 0) {

      cat("Potential mapping columns:\n")
      print(useful_cols)

      print(
        head(
          dat[, useful_cols, drop = FALSE],
          10
        )
      )
    }

  }, silent = TRUE)
}


# ------------------------------------------------------------
# 3. Search RDS objects for metadata containing TCGA IDs
# ------------------------------------------------------------

rds_files <- all_files[
  grepl(
    "\\.rds$",
    all_files,
    ignore.case = TRUE
  )
]

cat("\n============================================================\n")
cat(" SEARCHING RDS FILES\n")
cat("============================================================\n")

for (f in rds_files) {

  obj <- tryCatch(
    readRDS(f),
    error = function(e) NULL
  )

  if (is.null(obj)) next

  found_tcga <- FALSE
  found_individual <- FALSE

  # Data frame
  if (is.data.frame(obj)) {

    cols <- colnames(obj)

    found_tcga <- any(
      grepl(
        "TCGA-",
        unlist(
          lapply(
            obj,
            function(x) {
              if (is.character(x)) x else character(0)
            }
          )
        ),
        ignore.case = TRUE
      )
    )

    found_individual <- any(
      grepl(
        "individual",
        unlist(
          lapply(
            obj,
            function(x) {
              if (is.character(x)) x else character(0)
            }
          )
        ),
        ignore.case = TRUE
      )
    )

    if (found_tcga && found_individual) {

      cat("\n*** POSSIBLE MAPPING FOUND ***\n")
      cat("File:", f, "\n")
      cat("Class:", class(obj), "\n")
      cat("Dimensions:", dim(obj), "\n")
      print(head(obj, 10))

    }
  }
}


cat("\n============================================================\n")
cat(" SEARCH COMPLETE\n")
cat("============================================================\n")