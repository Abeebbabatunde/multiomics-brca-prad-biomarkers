# STEP 1.4 — Breast RPPA download and long-format file (GDC Proteome Profiling).
# Assembled from the code of the 19 June 2026 session: pair selection (before 12:31), download (12:31), long format + save (~13:00).
# Some statements were pasted twice; each appears once here, and console results are kept as comments. Nothing was added.
# Not in the source code: the creation of `clinical` is the GDCquery_clinic line from the earlier fragments (its columns submitter_id and bcr_patient_barcode were printed).
library(TCGAbiolinks)
query <- GDCquery(project = "TCGA-BRCA", data.category = "Proteome Profiling",
                  data.type = "Protein Expression Quantification")
results_df <- getResults(query)
nrow(results_df)                       # 919
table(results_df$sample_type)          # Metastatic 5, Primary Tumor 881, Solid Tissue Normal 33
clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")   # stage: early = Stage I, IA, IB, II, IIA, IIB (804 patients)
# ---- Select early-stage patients with both a Primary Tumor and a Solid Tissue Normal RPPA file ----
colnames(results_df)[grepl("submitter|case|patient", colnames(results_df), ignore.case = TRUE)]
# "cases" "submitter_id" "cases.submitter_id" "sample.submitter_id"
colnames(clinical)[grepl("submitter|case|patient", colnames(clinical), ignore.case = TRUE)]
# "submitter_id" "bcr_patient_barcode"

# Define early vs late stage
early_stages <- c("Stage I", "Stage IA", "Stage IB", "Stage II", "Stage IIA", "Stage IIB")
late_stages  <- c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IV")

clinical$stage_group <- ifelse(clinical$ajcc_pathologic_stage %in% early_stages, "Early",
                          ifelse(clinical$ajcc_pathologic_stage %in% late_stages, "Late", "Unknown"))

# Merge clinical stage info onto RPPA sample results
merged <- merge(results_df,
                 clinical[, c("submitter_id", "stage_group")],
                 by.x = "cases.submitter_id",
                 by.y = "submitter_id",
                 all.x = TRUE)

# Final filter: Normal samples + Early-stage tumors only
final_samples <- merged[merged$sample_type == "Solid Tissue Normal" |
                         (merged$sample_type == "Primary Tumor" & merged$stage_group == "Early"), ]

table(final_samples$sample_type, final_samples$stage_group, useNA = "always")
#                       Early Late <NA>
#   Primary Tumor         637    0    0
#   Solid Tissue Normal    30    3    0

final_samples$stage_group <- ifelse(final_samples$sample_type == "Solid Tissue Normal",
                                     "Normal", final_samples$stage_group)
table(final_samples$sample_type, final_samples$stage_group, useNA = "always")
#   Primary Tumor         637 Early, Solid Tissue Normal 33 Normal

normal_patients <- final_samples$cases.submitter_id[final_samples$sample_type == "Solid Tissue Normal"]
paired_tumor <- final_samples[final_samples$sample_type == "Primary Tumor" &
                               final_samples$cases.submitter_id %in% normal_patients, ]
nrow(paired_tumor)                     # 30

# Build final paired dataset: 30 patients, tumor + normal each
paired_patients <- intersect(normal_patients, paired_tumor$cases.submitter_id)

paired_normal <- final_samples[final_samples$sample_type == "Solid Tissue Normal" &
                                final_samples$cases.submitter_id %in% paired_patients, ]

paired_dataset <- rbind(paired_tumor, paired_normal)

table(paired_dataset$sample_type)
nrow(paired_dataset)  # should be 60 (30 tumor + 30 normal)

# Get the file/case identifiers for just our 60 paired samples
paired_barcodes <- paired_dataset$cases.submitter_id
paired_sample_types <- paired_dataset$sample_type
colnames(results_df)   # id, data_format, cases, access, file_name, submitter_id, data_category, type, platform, file_size, created_datetime, md5sum,
                       # updated_datetime, file_id, data_type, state, experimental_strategy, version, data_release, project, sample_type, is_ffpe,
                       # cases.submitter_id, sample.submitter_id

# Filter the query's internal results table to only the 60 paired samples
query$results[[1]] <- query$results[[1]][
  query$results[[1]]$cases.submitter_id %in% paired_dataset$cases.submitter_id,
]
nrow(query$results[[1]])               # 61
table(query$results[[1]]$sample_type)  # Metastatic 1, Primary Tumor 30, Solid Tissue Normal 30

query$results[[1]] <- query$results[[1]][query$results[[1]]$sample_type != "Metastatic", ]
table(query$results[[1]]$sample_type)
nrow(query$results[[1]])               # 60


GDCdownload(query)
rppa_data <- GDCprepare(query)         # 487 x 65 (AGID, lab_id, catalog_number, set_id, peptide_target + 60 samples)
metadata_cols <- c("AGID", "lab_id", "catalog_number", "set_id", "peptide_target")
sample_cols   <- setdiff(colnames(rppa_data), metadata_cols)       # 60
tumor_cols    <- sample_cols[grepl("-01A$|-01B$", sample_cols)]    # 30
normal_cols   <- sample_cols[grepl("-11A$|-11B$", sample_cols)]    # 30
protein_matrix_tumor  <- rppa_data[, c("peptide_target", tumor_cols)]
protein_matrix_normal <- rppa_data[, c("peptide_target", normal_cols)]
# Confirm no missing values that would break downstream stats (output not shown in the paste)
sum(is.na(protein_matrix_tumor))
sum(is.na(protein_matrix_normal))
dim(protein_matrix_tumor); dim(protein_matrix_normal)   # 487 x 31 each
# Confirm patient pairing actually lines up
tumor_patients  <- substr(tumor_cols, 1, 12)
normal_patients <- substr(normal_cols, 1, 12)
identical(sort(tumor_patients), sort(normal_patients))  # TRUE
# The first attempt, rbind(protein_matrix_tumor, protein_matrix_normal) after adding a sample_type column to each wide table, failed:
#   Error in match.names(clabs, names(xi)) : names do not match previous names
# Replaced by the pivot below (this is what produced the file). It needs library(tidyr), library(dplyr).
library(tidyr)
library(dplyr)

tumor_long <- protein_matrix_tumor %>%
  select(-sample_type) %>%  # remove the wide-format tag column added earlier by mistake
  pivot_longer(cols = -peptide_target, names_to = "sample_barcode", values_to = "value") %>%
  mutate(sample_type = "Primary Tumor")

normal_long <- protein_matrix_normal %>%
  select(-sample_type) %>%
  pivot_longer(cols = -peptide_target, names_to = "sample_barcode", values_to = "value") %>%
  mutate(sample_type = "Solid Tissue Normal")

brca_protein_raw <- bind_rows(tumor_long, normal_long)
dim(brca_protein_raw)  # 29220 x 4 (487 proteins x 60 samples)

# ---- Save the official Step 1.4 output ----
raw_dir <- "01_RawData/Breast"
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(brca_protein_raw, file.path(raw_dir, "brca_protein_raw.csv"), row.names = FALSE)
cat("Rows:", nrow(brca_protein_raw), "| Columns:", ncol(brca_protein_raw), "\n")

# Patient pairing reference (used by downstream steps)
pairing_info <- data.frame(
  patient_id = substr(unique(brca_protein_raw$sample_barcode), 1, 12),
  sample_barcode = unique(brca_protein_raw$sample_barcode),
  sample_type = brca_protein_raw$sample_type[match(unique(brca_protein_raw$sample_barcode), brca_protein_raw$sample_barcode)]
)
write.csv(pairing_info, file.path(raw_dir, "brca_protein_sample_pairing.csv"), row.names = FALSE)

# The same session also contained an audit for stray protein files and a check of the GDCdata/ cache size
# (suggested: unlink('GDCdata', recursive = TRUE) once the CSV is confirmed). Not part of the result chain.
