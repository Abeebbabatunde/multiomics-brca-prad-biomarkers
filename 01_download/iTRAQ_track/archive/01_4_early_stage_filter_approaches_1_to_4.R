# ================================================
# FILTER PROTEOMICS TO EARLY STAGE ONLY
# Run this after Step 1.4 download succeeds
# ================================================

cat("========================================\n")
cat("   FILTERING PROTEOMICS TO EARLY STAGE\n")
cat("========================================\n\n")

# -----------------------------------------------
# APPROACH 1: Match to your TCGA barcodes
# Best approach — keeps only your
# early stage patients from Step 1.1
# -----------------------------------------------
cat("APPROACH 1: Matching to TCGA barcodes...\n\n")

# Your early stage barcodes from Step 1.1
if(!exists("early_barcodes_brca")){
  clinical_brca <- read.csv(
    "01_RawData/Breast/clinical_brca.csv")
  early_barcodes_brca <- clinical_brca$submitter_id
}

cat("Your early-stage patients:",
    length(early_barcodes_brca), "\n\n")

# TCGA barcodes in proteomics are usually
# the first 12 characters of full barcode
# e.g. TCGA-A1-A0SB

# Shorten your barcodes to 12 characters
short_barcodes <- substr(
  early_barcodes_brca, 1, 12)

# Shorten protein column names
prot_short <- substr(
  colnames(brca_protein_raw), 1, 12)

# Find matching columns
matching_cols <- colnames(brca_protein_raw)[
  prot_short %in% short_barcodes
]

cat("Protein samples total    :",
    ncol(brca_protein_raw), "\n")
cat("Matching early-stage cols:",
    length(matching_cols), "\n\n")

if(length(matching_cols) >= 10){

  # Filter to early stage only
  brca_protein_early <- brca_protein_raw[,
    matching_cols]

  cat("✅ APPROACH 1 SUCCESS!\n")
  cat("  Early-stage protein samples:",
      ncol(brca_protein_early), "\n\n")

} else {

  cat("⚠️  Few barcode matches found\n")
  cat("    Protein sample names may use\n")
  cat("    different format than TCGA\n")
  cat("    Trying Approach 2...\n\n")

  # -----------------------------------------------
  # APPROACH 2: Filter using clinical stage
  # in the proteomics metadata itself
  # -----------------------------------------------
  cat("APPROACH 2: Using proteomics metadata...\n\n")

  if(exists("brca_prot_meta")){

    cat("Sample types in proteomics data:\n")
    if("disease state:ch1" %in%
       colnames(brca_prot_meta)){
      print(table(brca_prot_meta$
                  `disease state:ch1`))
    }
    if("characteristics_ch1" %in%
       colnames(brca_prot_meta)){
      print(table(brca_prot_meta$
                  characteristics_ch1))
    }
    if("ajcc_pathologic_stage" %in%
       colnames(brca_prot_meta)){
      print(table(brca_prot_meta$
                  ajcc_pathologic_stage))
    }

    cat("\nAll metadata columns:\n")
    print(colnames(brca_prot_meta))

  } else {
    cat("⚠️  No metadata available\n")
    cat("    Using Approach 3...\n\n")
  }
}

# -----------------------------------------------
# APPROACH 3: If proteomics has stage info
# Filter directly using stage column
# -----------------------------------------------

filter_by_stage <- function(meta_df,
                              data_mat){

  # Find stage column
  stage_cols <- colnames(meta_df)[
    grepl("stage|Stage|STAGE|ajcc|tumor",
          colnames(meta_df),
          ignore.case=TRUE)
  ]

  cat("Stage-related columns found:\n")
  print(stage_cols)

  if(length(stage_cols) > 0){

    stage_col <- stage_cols[1]
    cat("\nUsing column:", stage_col, "\n")
    cat("Stage values:\n")
    print(table(meta_df[[stage_col]]))

    # Keep Stage I and II
    early_idx <- grepl(
      "Stage I|Stage II|^I$|^II$|^IA|^IB|^IIA|^IIB",
      meta_df[[stage_col]],
      ignore.case=TRUE
    )

    early_samples <- rownames(
      meta_df)[early_idx]

    common <- intersect(
      early_samples,
      colnames(data_mat)
    )

    if(length(common) >= 10){
      cat("\n✅ Early-stage samples found:",
          length(common), "\n")
      return(data_mat[, common])
    } else {
      cat("\n⚠️  Too few early-stage samples:",
          length(common), "\n")
      return(NULL)
    }
  } else {
    cat("No stage column found\n")
    return(NULL)
  }
}

# Try stage filtering if approach 1 failed
if(!exists("brca_protein_early") &&
   exists("brca_prot_meta")){

  result <- filter_by_stage(
    brca_prot_meta,
    brca_protein_raw
  )

  if(!is.null(result)){
    brca_protein_early <- result
    cat("✅ APPROACH 2/3 SUCCESS!\n")
    cat("  Early-stage samples:",
        ncol(brca_protein_early), "\n\n")
  }
}

# -----------------------------------------------
# APPROACH 4: If no stage info available
# Use ALL samples but label them clearly
# This is still scientifically valid
# -----------------------------------------------
if(!exists("brca_protein_early")){

  cat("APPROACH 4: Using all available samples\n")
  cat("with clear documentation...\n\n")

  brca_protein_early <- brca_protein_raw

  cat("⚠️  NOTE FOR YOUR THESIS:\n")
  cat("   The proteomics dataset does not\n")
  cat("   have explicit stage labels.\n")
  cat("   You will document this as a\n")
  cat("   limitation and note that RNA-seq\n")
  cat("   and methylation are confirmed\n")
  cat("   early-stage (Stage I and II)\n\n")
}

# -----------------------------------------------
# SHOW WHAT YOU HAVE
# -----------------------------------------------
cat("========================================\n")
cat("   YOUR PROTEOMICS DATASET SUMMARY\n")
cat("========================================\n\n")

cat(sprintf("  Total proteins   : %d\n",
            nrow(brca_protein_early)))
cat(sprintf("  Total samples    : %d\n",
            ncol(brca_protein_early)))

# Check for cancer vs normal labels
if(exists("brca_prot_meta")){
  type_cols <- colnames(brca_prot_meta)[
    grepl("type|status|normal|tumor|cancer",
          colnames(brca_prot_meta),
          ignore.case=TRUE)
  ]
  if(length(type_cols) > 0){
    cat(sprintf("  Sample types     :\n"))
    print(table(brca_prot_meta[[type_cols[1]]]))
  }
}