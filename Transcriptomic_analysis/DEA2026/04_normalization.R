#########################################################################################################
# Normalized (batch-corrected) expression matrix + metadata, aligned to the same
# samples/genes used in 02_rnaseq_2026.R's DESeq2 run.
#########################################################################################################

start_time <- Sys.time()
cat("Start time is: ", format(start_time), "\n")

library("limma")
library("DESeq2")
library("dplyr")
library("edgeR")

################################################################################
################################## Read data ##################################
################################################################################

# Same inputs as 02_rnaseq_2026.R -- guarantees this normalization matches the
# exact sample/gene set that went into the DESeq2 model.
countData <- readRDS(file = "./data/rawCounts_filtered.rds")
metadata  <- readRDS(file = "./data/confoundingDataFrame2_full.rds")


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
######################## TMM normalization + logCPM ###########################
################################################################################

dge <- DGEList(counts = countData)

normFactors <- calcNormFactors(countData, method = "TMM", refColumn = 1,
                               logratioTrim = 0.3, sumTrim = 0.05,
                               doWeighting = TRUE, Acutoff = -1e+10)
normFactors <- normFactors * (colSums(countData) / mean(colSums(countData)))
names(normFactors) <- colnames(countData)

logCPM <- cpm(dge, log = TRUE, prior.count = 1)

################################################################################
######################## Batch correction (nuisance vars only) ################
################################################################################

# Column renames vs. old normalization_junio.R:
#   condition -> case_control (kept OUT of batch correction -- see below)
#   neutPer   -> Neutrophils
# case_control and visit_month are deliberately excluded from selectedVars:
# they are variables of biological interest, not nuisance covariates, and
# should not be regressed out of the continuous expression matrix.

selectedVars <- c("Neutrophils", "RIN_Value", "Plate", "Concentration",
                  "sex", "age_at_baseline", "race")

for (v in c("Submitted_Volume__ul_", "Concentration")) {
  if (!v %in% colnames(metadata)) {
    warning(v, " not found in metadata -- omitted from batch correction. ",
            "Confirm this column is genuinely absent, not just renamed.")
    selectedVars <- setdiff(selectedVars, v)
  } else if (!v %in% selectedVars) {
    selectedVars <- c(selectedVars, v)
  }
}

na_check <- sapply(metadata[, selectedVars], function(x) sum(is.na(x)))
print(na_check)
if (sum(na_check) > 0) {
  stop("NA values present in batch-correction covariates -- model.matrix() ",
       "would silently drop these samples. Fix upstream before proceeding.")
}

matrix2Correct <- metadata[, selectedVars]
designCovariates <- model.matrix(~., data = matrix2Correct)[, -1]

cat("Is batch-correction design fullrank?\n")
print(is.fullrank(model.matrix(~., data = matrix2Correct)))

normExpr <- removeBatchEffect(logCPM, covariates = designCovariates)

################################################################################
################################## Save outputs ###############################
################################################################################


saveRDS(metadata, file = "normalized_metadata.rds")
saveRDS(normExpr, file = "normalized_expression.rds")

cat("Normalization done. Samples:", ncol(normExpr), "| Genes:", nrow(normExpr), "\n")

end_time <- Sys.time()
cat("Time taken for analysis: ", format(end_time - start_time), "\n")