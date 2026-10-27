################################################################################
# Differential expression (DESeq2) -- PD cases vs controls, whole-blood RNA-seq
# Input: checkpoints saved by 01_data_processing_rna_2026.R
################################################################################

start_time <- Sys.time()
cat("Start time is: ", format(start_time), "\n")

library(DESeq2)
library(dplyr)
library(edgeR)
library(BiocParallel)

register(MulticoreParam(workers = 8))   # adjust to the machine

setwd("~/Documents/Omics_Integration/Transcriptomic_analysis/DEA2026/")
dir.create("results/Results_rnaseq", recursive = TRUE, showWarnings = FALSE)

################################################################################
## 1. Read and align data
################################################################################

countData <- readRDS("./data/rawCounts_filtered.rds")
metadata  <- droplevels(readRDS("./data/confoundingDataFrame2_full.rds"))
# Gene filtering and MT/ribo/globin exclusion were done in preprocessing -- do NOT redo here.

common_samples <- intersect(colnames(countData), rownames(metadata))
countData <- countData[, common_samples, drop = FALSE]
metadata  <- metadata[common_samples, , drop = FALSE]

stopifnot(identical(colnames(countData), rownames(metadata)),
          ncol(countData) == 5731,                 # expected N after Plate x race filter
          !anyDuplicated(rownames(countData)))     # PAR duplicates must not reach the model
cat("Samples:", ncol(countData), "| Genes:", nrow(countData), "\n")

################################################################################
## 2. Design
################################################################################

design_terms <- c("case_control", "visit_month", "sex", "race", "age_at_baseline",
                  "Neutrophils", "Plate", "RIN_Value", "Normalization_Volume__30ng_ul_")
stopifnot(all(design_terms %in% colnames(metadata)))

na_check <- sapply(metadata[, design_terms], function(x) sum(is.na(x)))
print(na_check)
stopifnot(sum(na_check) == 0)

# Reference = control, so log2FC > 0 means higher in cases
metadata$case_control <- relevel(factor(metadata$case_control), ref = "control")
metadata$sex   <- factor(metadata$sex)
metadata$race  <- factor(metadata$race)
metadata$Plate <- factor(metadata$Plate)

# Scale continuous covariates (better GLM convergence; does not change the diagnosis estimate)
cont_vars <- c("age_at_baseline", "Neutrophils", "RIN_Value", "Normalization_Volume__30ng_ul_")
metadata[cont_vars] <- lapply(metadata[cont_vars], function(x) as.numeric(scale(x)))

for (v in design_terms) {
  n_levels <- if (is.factor(metadata[[v]])) nlevels(metadata[[v]]) else length(unique(metadata[[v]]))
  cat(v, ": ", n_levels, " levels/unique values\n", sep = "")
}

design <- model.matrix(as.formula(paste("~", paste(design_terms, collapse = " + "))), metadata)
stopifnot(nrow(design) == ncol(countData),   # catches silent NA-row dropping
          limma::is.fullrank(design))
cat("Design:", nrow(design), "samples x", ncol(design), "coefficients\n")

################################################################################
## 3. TMM-derived size factors (same settings as previous run)
################################################################################

normFactors <- calcNormFactors(countData, method = "TMM", refColumn = 1,
                               logratioTrim = 0.3, sumTrim = 0.05,
                               doWeighting = TRUE, Acutoff = -1e+10)
sizeFactors_tmm <- normFactors * (colSums(countData) / mean(colSums(countData)))
names(sizeFactors_tmm) <- colnames(countData)

################################################################################
## 4. DESeq2
################################################################################

dds <- DESeqDataSetFromMatrix(countData = countData, colData = metadata, design = design)
sizeFactors(dds) <- sizeFactors_tmm

cat("Running DESeq()\n")
dds <- DESeq(dds, parallel = TRUE)
saveRDS(dds, "results/Results_rnaseq/dds_2026.rds")

coef_dx <- grep("^case_control", resultsNames(dds), value = TRUE)
stopifnot(length(coef_dx) == 1)
cat("Diagnosis coefficient:", coef_dx, "(log2FC > 0 = higher in cases)\n")

res <- as.data.frame(results(dds, name = coef_dx))
res <- res[order(res$padj), ]
saveRDS(res, "results/Results_rnaseq/res_deseq_2026.rds")
write.table(res, "results/Results_rnaseq/genesDiffExpressedParkinson.tsv",
            sep = "\t", quote = FALSE, col.names = NA)

################################################################################
## 5. Summary
################################################################################

degs <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 0.15, ]
cat("DEGs (padj < 0.05, |log2FC| > 0.15):", nrow(degs),
    "| higher in cases:", sum(degs$log2FoldChange > 0),
    "| lower in cases:", sum(degs$log2FoldChange < 0), "\n")
print(head(degs[, c("baseMean", "log2FoldChange", "padj")], 30))

writeLines(capture.output(sessionInfo()), "results/Results_rnaseq/sessionInfo_deseq2_2026.txt")

end_time <- Sys.time()
cat("Time taken for analysis: ", format(end_time - start_time), "\n")