SUPERSEDED -- do not run for the thesis results.
These scripts built the earlier tumour-only integration cohorts (31 BRCA three-omics patients; 187 PRAD two-omics patients).
The reported analysis uses the paired Normal-vs-Early design in ../../05_0a..05_6 (BRCA 51/16 pairs, PRAD 17 pairs).
Facts these logs establish:
 - 00.5 manifest: 991 patients x 14 columns (BRCA 803, PRAD 188). Attrition report: BRCA epi 570, txn 801, prot 76, complete 31; PRAD epi 188, txn 187, prot 0, complete 187.
 - The BRCA three-omics 31 were tumour-only; 3 had no 01 methylation sample (TCGA-A7-A0CE, TCGA-BH-A0BV = code 11; TCGA-BH-A18V = code 06).
 - PRAD: 190 primary-tumour (01) columns in both RNA and methylation for the 187 patients (3 patients have 01A+01B: TCGA-HC-7740, TCGA-HC-8258, TCGA-HC-8261).
 - PRAD proteomics (36 samples, "individual N, neoplasm / normal tissue adjacent to tumor") has no TCGA IDs; analysed separately (the match-count output was truncated in the log).
Missing from this folder (pasted into chat, not saved as files): the 5.0C / 5.0E searches, the 5.2A and 5.3F/5.3G script text (the console copies are here).
