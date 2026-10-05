################################################################################
# Standalone diagnostic: swamp / PCA on raw Olink proteomics data (CSF + PLA)
# Independent of the limma DE scripts -- reads raw data itself, only writes
# to Results/Bias_correction_raw/ so it can never overwrite DE outputs or
# any diagnostic files the limma scripts might later produce.
################################################################################

library(dplyr)
library(swamp)
library(ggplot2)

setwd("~/Documents/Parkinson/proteomic/")
dir.create("~/Documents/Omics_Integration/Proteomic_analysis/Results_full_aligned/Bias_correction_raw",
           recursive = TRUE, showWarnings = FALSE)
out_dir <- "~/Documents/Omics_Integration/Proteomic_analysis/Results_full_aligned/Bias_correction_raw"

# ------------------------------------------------------------------------- #
# Reusable loader: reads one matrix type (CSF or PLA) from its Olink folder
# ------------------------------------------------------------------------- #
load_olink_matrix <- function(matrix_dir_tag) {
  base <- paste0("Targeted olink/proteomics-", matrix_dir_tag, "-PPEA-D03/")

  samples_full <- read.csv(file = paste0(base, "releases_2023_v4release_1027_proteomics-",
                                          matrix_dir_tag, "-PPEA-D03_proteomics_",
                                          matrix_dir_tag, "-PPEA-D03_samples.csv"))

  panel_dir <- paste0(base, "protein-expression/releases_2023_v4release_1027_proteomics-",
                       matrix_dir_tag, "-PPEA-D03_olink-explore_protein-expression_",
                       matrix_dir_tag, "-PPEA-D03_matrix_")

  cardiometa   <- read.csv(file = paste0(panel_dir, "cardiometabolic.csv"))
  inflammatory <- read.csv(file = paste0(panel_dir, "inflammation.csv"))
  neurology    <- read.csv(file = paste0(panel_dir, "neurology.csv"))
  oncology     <- read.csv(file = paste0(panel_dir, "oncology.csv"))

  proteomic_raw <- rbind(cardiometa, inflammatory, neurology, oncology)
  rownames(proteomic_raw) <- proteomic_raw[[1]]
  proteomic_raw <- proteomic_raw[, -1]
  stopifnot(!anyDuplicated(rownames(proteomic_raw)))

  metadata <- readRDS(file = "~/Documents/Parkinson/proteomic/metadata.rds")

  samples <- merge(samples_full, metadata[, c(1, 11, 13, 14)])
  colnames(samples)[colnames(samples) == "case_control_other_at_baseline"] <- "condition"
  samples <- samples %>% dplyr::select(participant_id, sample_id, condition, sex, age_at_baseline, visit_month)
  samples$condition <- factor(samples$condition)
  samples$sample_id <- gsub("-", ".", samples$sample_id)
  rownames(samples) <- samples$sample_id

  full_normalized_counts <- as.matrix(proteomic_raw)

  list(matrix = full_normalized_counts, samples = samples)
}

# ------------------------------------------------------------------------- #
# Diagnostic: swamp::prince() + PC-pair scatter plots, PNG + PDF
# ------------------------------------------------------------------------- #
run_swamp_diagnostic <- function(full_normalized_counts, samples, tag) {
  stopifnot(all(samples$sample_id %in% colnames(full_normalized_counts)))
  diag_meta <- samples[colnames(full_normalized_counts), c("condition", "visit_month", "sex", "age_at_baseline")]
  diag_meta[] <- lapply(diag_meta, function(x) if (is.character(x)) factor(x) else x)
  stopifnot(all(rownames(diag_meta) == colnames(full_normalized_counts)))

  # Olink NPX can have NAs below LOD -- prcomp/prince need complete cases.
  na_per_protein <- rowSums(is.na(full_normalized_counts))
  diag_matrix <- full_normalized_counts[na_per_protein == 0, ]
  cat("[", tag, "] Proteins kept for diagnostic (complete cases):", nrow(diag_matrix),
      "of", nrow(full_normalized_counts), "\n")
  stopifnot(nrow(diag_matrix) > 0)

  # Shared helper to save a prince() result as swamp PNG + PDF
  save_swamp <- function(resConf, label) {
    saveRDS(resConf, file.path(out_dir, paste0("resConf_proteomics_", tag, "_", label, ".rds")))

    png(file.path(out_dir, paste0("swamp_proteomics_", tag, "_", label, ".png")), width = 960, height = 960)
    prince.plot(resConf, margins = c(5, 15))
    dev.off()

    pdf(file.path(out_dir, paste0("swamp_proteomics_", tag, "_", label, ".pdf")), width = 10, height = 10)
    prince.plot(resConf, margins = c(5, 15))
    dev.off()
  }

  # -------------------------------------------------------------------- #
  # BEFORE correction
  # -------------------------------------------------------------------- #
  resConf_before <- prince(diag_matrix, diag_meta, top = 10)
  save_swamp(resConf_before, "before")

  # -------------------------------------------------------------------- #
  # Diagnostic-only correction: adjust for known confounders EXCLUDING
  # `condition` (the variable of interest). Never fed into any DE model --
  # purely to visualize whether these confounders explain PC variance.
  # -------------------------------------------------------------------- #
  correctVars <- c("visit_month", "sex", "age_at_baseline")
  correctVars <- correctVars[correctVars %in% colnames(diag_meta)]
  na_check <- sapply(diag_meta[, correctVars], function(x) sum(is.na(x)))
  stopifnot(sum(na_check) == 0)

  designCovariates <- model.matrix(~., data = diag_meta[, correctVars, drop = FALSE])[, -1]
  diag_matrix_corrected <- limma::removeBatchEffect(diag_matrix, covariates = designCovariates)

  # -------------------------------------------------------------------- #
  # AFTER correction
  # -------------------------------------------------------------------- #
  resConf_after <- prince(as.matrix(diag_matrix_corrected), diag_meta, top = 10)
  save_swamp(resConf_after, "after")

  cat("[", tag, "] swamp before/after plots saved to", out_dir, "\n")
}

# ------------------------------------------------------------------------- #
# Run for both matrices
# ------------------------------------------------------------------------- #
csf_data <- load_olink_matrix("CSF")
run_swamp_diagnostic(csf_data$matrix, csf_data$samples, "CSF")

pla_data <- load_olink_matrix("PLA")
run_swamp_diagnostic(pla_data$matrix, pla_data$samples, "PLA")