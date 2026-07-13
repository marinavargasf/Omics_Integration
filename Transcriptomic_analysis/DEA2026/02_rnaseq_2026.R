################################################################################
################################################################################

start_time <- Sys.time()
cat("Start time is: ", format(start_time), "\n")

library("DESeq2")
library("dplyr")
library("edgeR")

################################################################################
################################## Read data ##################################
################################################################################
#setwd("~/Documents/Omics_Integration/Transcriptomic_analysis/DEA2026/")

countData <- readRDS(file = "./data/rawCounts_filtered.rds")
metadata   <- readRDS(file = "./data/confoundingDataFrame2_full.rds")

# Gene filtering, MT/ribo/globin exclusion, and Ensembl mapping were already
# done in preprocessing (with Ensembl pinned to v100) -- do NOT redo them here,
# redoing with an unpinned useMart() call is what caused the DEG list drift.

################################################################################
###################### Align countData <-> metadata ##########################
################################################################################

common_samples <- intersect(colnames(countData), rownames(metadata))
stopifnot(length(common_samples) > 0)
countData <- countData[, common_samples, drop = FALSE]
metadata  <- metadata[common_samples, , drop = FALSE]
countData <- countData[, rownames(metadata)]

stopifnot(all(rownames(metadata) == colnames(countData)))
cat("Samples aligned:", ncol(countData), "\n")

################################################################################
################# Differential expression on RNA-seq data #####################
################################################################################

normFactors <- calcNormFactors(countData, method = "TMM", refColumn = 1,
                               logratioTrim = 0.3, sumTrim = 0.05,
                               doWeighting = TRUE, Acutoff = -1e+10)
normFactors <- normFactors * (colSums(countData) / mean(colSums(countData)))
names(normFactors) <- colnames(countData)

cat("Normalization factors calculated\n")

# NOTE column renames vs. the old script:
#   condition -> case_control   (check resultsNames(dds) below for the actual
#                                 coefficient name -- don't assume it matches
#                                 the old "conditionControl")
#   neutPer   -> Neutrophils
# X_260_280_Ratio dropped: preprocessing removes all "*Ratio*" QC columns.
# Submitted_Volume__ul_ kept only if present in metadata -- checked below.

design_terms <- c("case_control", "visit_month", "sex", "age_at_baseline",
                  "race", "Neutrophils", "Plate", "RIN_Value")
if ("Submitted_Volume__ul_" %in% colnames(metadata)) {
  design_terms <- c(design_terms, "Submitted_Volume__ul_")
} else {
  warning("Submitted_Volume__ul_ not found in metadata -- omitted from design. ",
          "Confirm this column is genuinely gone from the QC table, not just misnamed.")
}

na_check <- sapply(metadata[, design_terms], function(x) sum(is.na(x)))
print(na_check)
if (sum(na_check) > 0) {
  stop("NA values present in design covariates -- model.matrix() would silently ",
       "drop these samples. Fix upstream before proceeding.")
}

design_formula <- as.formula(paste("~", paste(design_terms, collapse = " + ")))
design <- model.matrix(design_formula, metadata)

cat("Is design fullrank? \n")
print(is.fullrank(design))

cat("Does metadata and countData match? must be TRUE\n")
print(all(rownames(metadata) == colnames(countData)))

cat("Metadata dim\n")
print(dim(metadata))
cat("countData dim\n")
print(dim(countData))
cat("design ncol\n")
print(ncol(design))

stopifnot(nrow(design) == ncol(countData))  # catches silent NA-row-dropping

################################################################################

ddsFromMatrix <- DESeqDataSetFromMatrix(countData = countData,
                                        colData = metadata,
                                        design = design)

cat("Assigning TMM-derived size factors\n")
sizeFactors(ddsFromMatrix) <- normFactors

cat("Running DESeq()\n")
dds <- DESeq(ddsFromMatrix, parallel = TRUE)

cat("DESeq done!\n")
cat("Available result coefficients:\n")
print(resultsNames(dds))

save.image("Results_rnaseq/deseq2_output.RData")

# Uncomment once you've confirmed the correct coefficient name from resultsNames(dds):
# res <- as.data.frame(results(dds, name = "case_control_control_vs_case"))
# saveRDS(res, file = "Results_2026/res_deseq.rds")

end_time <- Sys.time()
cat("Time taken for analysis: ", format(end_time - start_time), "\n")