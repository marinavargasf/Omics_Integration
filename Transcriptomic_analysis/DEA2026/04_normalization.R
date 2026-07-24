#########################################################################################################
# Normalized (batch-corrected) expression matrix + metadata, aligned to the same
# samples/genes used in 02_rnaseq_2026.R's DESeq2 run.
#########################################################################################################
setwd("~/Documents/Omics_Integration/Transcriptomic_analysis/DEA2026/")
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


################################################################################
##################### Translate to gene symbol (HGNC) #########################
################################################################################

#normExpr <- readRDS("~/Documents/Omics_Integration/Transcriptomic_analysis/DEA2026/results/Results_rnaseq/normalized_expression.rds")

ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

ensembl_ids <- rownames(normExpr)
gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

annotations <- data.frame(ensembl_gene_id = ensembl_ids, gene_name = NA, stringsAsFactors = FALSE)
match_idx <- match(annotations$ensembl_gene_id, gene_mapping$ensembl_gene_id)
annotations$gene_name <- gene_mapping$hgnc_symbol[match_idx]

stopifnot(nrow(annotations) == nrow(normExpr))

normExpr_df <- as.data.frame(normExpr)
normExpr_df$gene_name <- annotations$gene_name
missing_symbol <- is.na(normExpr_df$gene_name) | normExpr_df$gene_name == ""
normExpr_df$gene_name[missing_symbol] <- rownames(normExpr_df)[missing_symbol]

# duplicate gene symbols: keep the row with highest mean expression across samples
# (alternative: highest variance -- more relevant if this feeds a correlation/
# network step where variability matters more than absolute level; swap the
# ranking metric below if that's a better fit for MariNET)
expr_cols <- setdiff(colnames(normExpr_df), "gene_name")
normExpr_df$mean_expr <- rowMeans(normExpr_df[, expr_cols])

data_ordered <- normExpr_df[order(normExpr_df$mean_expr, decreasing = TRUE), ]
data_ordered <- data_ordered[!duplicated(data_ordered$gene_name), ]
data_ordered <- data_ordered[!is.na(data_ordered$gene_name), ]

rownames(data_ordered) <- data_ordered$gene_name
data_ordered <- data_ordered[, expr_cols]

normExpr_genesymbol <- as.matrix(data_ordered)

saveRDS(normExpr_genesymbol, file = "./results/Results_rnaseq/normalized_expression_genesymbol.rds")

n_unmapped <- sum(missing_symbol)
cat("Gene symbol translation done. Genes:", nrow(normExpr_genesymbol),
    "(from", nrow(normExpr), "Ensembl IDs;", n_unmapped, "had no HGNC symbol and kept their Ensembl ID;",
    "duplicate symbols were collapsed to the highest-mean-expression row)\n")

end_time <- Sys.time()
cat("Time taken for analysis: ", format(end_time - start_time), "\n")
