################################################################################
# 00a_setup.R
# Run once before the pipeline (from the repository root / .Rproj).
#   1. Checks that the required R packages are installed (and MariNET >= 2.0.0)
#   2. Creates the output folder tree under results/
#   3. Checks that the external inputs (AMP PD / PPMI data and resources/)
#      are in place and reports anything missing
# It does not modify any data.
################################################################################

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
setwd(here::here())

## 1. Packages ----------------------------------------------------------------

cran_pkgs <- c("here", "dplyr", "tidyr", "purrr", "tibble", "readr", "data.table",
               "ggplot2", "ggrepel", "scales", "pheatmap", "rstatix", "ggpubr",
               "lme4", "qgraph", "ggeffects", "MuMIn", "swamp", "msigdbr", "MariNET",
               "plyr", "table1", "psych", "bootnet")
bioc_pkgs <- c("DESeq2", "edgeR", "limma", "biomaRt", "Biobase", "BiocParallel",
               "fgsea", "pathMED")

missing_pkgs <- c(cran_pkgs, bioc_pkgs)[!vapply(c(cran_pkgs, bioc_pkgs),
                                                requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall them with:\n  install.packages('BiocManager')\n",
       "  BiocManager::install(c(", paste0('"', missing_pkgs, '"', collapse = ", "), "))")
}
stopifnot("MariNET >= 2.0.0 is required (lmm_integration_analysis)" =
            packageVersion("MariNET") >= "2.0.0")
cat("All required packages are installed.\n")

## 2. Output folders ----------------------------------------------------------

out_dirs <- c("results/clinical",
              "results/transcriptomics/data",
              "results/transcriptomics/bias_correction",
              "results/transcriptomics/deseq2",
              "results/transcriptomics/marinet",
              "results/pathways",
              "results/proteomics/bias_correction",
              "results/proteomics/DDC_UPDRS4_followup")
for (d in out_dirs) dir.create(here::here(d), recursive = TRUE, showWarnings = FALSE)
cat("Output folders created under results/.\n")

## 3. External inputs ---------------------------------------------------------

amp  <- "data/raw/AMP_PD_Data_v3_2022"
olink <- "data/raw/proteomics/Targeted olink"
rel  <- "releases_2023_v4release_1027_"

required <- c(
  # AMP PD v4.0 - clinical tables (clinical_processing/clinical_data.R reads the whole folder)
  file.path(amp, "clinical_metadata_v4_2023", paste0(rel, "amp_pd_participants.csv")),
  file.path(amp, "clinical_metadata_v4_2023/clinical", paste0(rel, "clinical_Demographics.csv")),
  file.path(amp, "clinical_metadata_v4_2023/clinical", paste0(rel, "clinical_PD_Medical_History.csv")),
  file.path(amp, "clinical_metadata_v4_2023", paste0(rel, "amp_pd_case_control.csv")),
  # AMP PD v4.0 - whole-blood RNA-seq
  file.path(amp, "transcriptomics_v4/reports", paste0(rel, "rnaseq-WB-RWTS_sequencing_metrics_rna_quality_metrics.csv")),
  file.path(amp, "transcriptomics_v4/reports", paste0(rel, "rnaseq-WB-RWTS_reports_multiqc_multiqc_data_multiqc_star.txt")),
  file.path(amp, "transcriptomics_v4/subread_feature-counts", paste0(rel, "rnaseq-WB-RWTS_subread_feature-counts_matrix.featureCounts.tsv")),
  # PPMI hematology (neutrophil counts)
  "data/raw/PPMI/Blood_Chemistry___Hematology_10Mar2025.csv",
  # Fixed reference files (see resources/README.md)
  "resources/clinical_layout.rds",
  "resources/uniprot_to_hgnc_mapping_full.rds"
)
for (tissue in c("CSF", "PLA")) {
  base <- file.path(olink, paste0("proteomics-", tissue, "-PPEA-D03"))
  required <- c(required,
    file.path(base, paste0(rel, "proteomics-", tissue, "-PPEA-D03_proteomics_", tissue, "-PPEA-D03_samples.csv")),
    file.path(base, paste0(rel, "proteomics-", tissue, "-PPEA-D03_ColumnHeaderGlossary.csv")),
    file.path(base, "protein-expression",
              paste0(rel, "proteomics-", tissue, "-PPEA-D03_olink-explore_protein-expression_",
                     tissue, "-PPEA-D03_matrix_",
                     c("cardiometabolic", "inflammation", "neurology", "oncology"), ".csv")))
}

# Optional: cached Ensembl annotations. If absent they are rebuilt from the
# current Ensembl release, which may differ slightly from the published run.
optional <- c("resources/gene_map_v100_full.rds",
              "resources/exclude_genes_ENSEMBL_v100.rds")

ok_req <- file.exists(here::here(required))
ok_opt <- file.exists(here::here(optional))

if (any(!ok_opt)) {
  message("Optional cached annotation(s) not found (will be rebuilt from Ensembl):\n  ",
          paste(optional[!ok_opt], collapse = "\n  "))
}
if (any(!ok_req)) {
  stop(sum(!ok_req), " required input(s) not found:\n  ",
       paste(required[!ok_req], collapse = "\n  "),
       "\nSee README.md, section 'Input data'.")
}
cat("All required inputs found. You can run 00b and then scripts 01-16 in order.\n")
