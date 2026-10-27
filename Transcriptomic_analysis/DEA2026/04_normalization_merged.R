#########################################################################################################
# Normalized (batch-corrected) expression matrix + metadata, aligned to the same
# samples/genes used in 02_deseq2_2026.R's DESeq2 run. For MariNET network input.
#########################################################################################################
setwd("~/Documents/Omics_Integration/Transcriptomic_analysis/DEA2026/results")
start_time <- Sys.time()
cat("Start time is: ", format(start_time), "\n")

library("limma")
library("DESeq2")
library("dplyr")
library("edgeR")

################################################################################
################################## Read data ##################################
################################################################################

# Same inputs as 02_deseq2_2026.R -- already gene-filtered, QC'd, rare-combo
# filtered, droplevels()'d. Do NOT redo any of that here.
# NOTE: working dir is results/, so data/ is one level up.
countData <- readRDS(file = "../data/rawCounts_filtered.rds")
metadata  <- droplevels(readRDS(file = "../data/confoundingDataFrame2_full.rds"))

common_samples <- intersect(colnames(countData), rownames(metadata))
stopifnot(length(common_samples) > 0)
countData <- countData[, common_samples, drop = FALSE]
metadata  <- metadata[common_samples, , drop = FALSE]
stopifnot(identical(rownames(metadata), colnames(countData)),
          ncol(countData) == 5731)                 # same N as the DESeq2 run
cat("Samples aligned:", ncol(countData), "\n")

################################################################################
######################## TMM normalization + logCPM ###########################
################################################################################

dge <- DGEList(counts = countData)
dge <- calcNormFactors(dge, method = "TMM")
logCPM <- cpm(dge, log = TRUE, prior.count = 1)

################################################################################
######################## Batch correction (nuisance vars only) ################
################################################################################

# Same covariates as the DESeq2 design (02_deseq2_2026.R) minus case_control and
# visit_month: those are biological signal for MariNET and are PRESERVED via `design`.
selectedVars <- c("Neutrophils", "RIN_Value", "Plate", "sex", "age_at_baseline", "race",
                  "Normalization_Volume__30ng_ul_")
stopifnot(all(selectedVars %in% colnames(metadata)))   # fail loudly instead of silently dropping
cat("Batch-correction covariates used:", paste(selectedVars, collapse = ", "), "\n")

na_check <- sapply(metadata[, c(selectedVars, "case_control", "visit_month")], function(x) sum(is.na(x)))
print(na_check)
stopifnot(sum(na_check) == 0)

matrix2Correct <- droplevels(metadata[, selectedVars])
designCovariates <- model.matrix(~., data = matrix2Correct)[, -1]
stopifnot(is.fullrank(model.matrix(~., data = matrix2Correct)))

design2preserve <- model.matrix(~ case_control + visit_month, data = metadata)

normExpr <- removeBatchEffect(logCPM, covariates = designCovariates, design = design2preserve)

################################################################################
##################### Annotate with gene symbols (cached mapping) ############
################################################################################

gene_map_path <- "../data/gene_map_v100_full.rds"
stopifnot(file.exists(gene_map_path))  # built in 03_annotate_deseq2_2026.R -- reuse it
gene_mapping <- readRDS(gene_map_path)
cat("Loaded cached gene mapping (", nrow(gene_mapping), "rows)\n")

ensembl_ids <- rownames(normExpr)
match_idx <- match(ensembl_ids, gene_mapping$ensembl_gene_id)
gene_name <- gene_mapping$hgnc_symbol[match_idx]
missing_symbol <- is.na(gene_name) | gene_name == ""
gene_name[missing_symbol] <- ensembl_ids[missing_symbol]

normExpr_df <- as.data.frame(normExpr)
normExpr_df$gene_name <- gene_name
normExpr_df$mean_expr <- rowMeans(normExpr)

data_ordered <- normExpr_df[order(normExpr_df$mean_expr, decreasing = TRUE), ]
data_ordered <- data_ordered[!duplicated(data_ordered$gene_name), ]
rownames(data_ordered) <- data_ordered$gene_name
normExpr_genesymbol <- as.matrix(data_ordered[, setdiff(colnames(normExpr_df), c("gene_name", "mean_expr"))])

################################################################################
################################## Save outputs ###############################
################################################################################

dir.create("Results_rnaseq", recursive = TRUE, showWarnings = FALSE)
saveRDS(metadata, file = "Results_rnaseq/normalized_metadata.rds")
saveRDS(normExpr, file = "Results_rnaseq/normalized_expression.rds")
saveRDS(normExpr_genesymbol, file = "Results_rnaseq/normalized_expression_genesymbol.rds")

cat("Gene symbol translation done. Genes:", nrow(normExpr_genesymbol),
    "(from", nrow(normExpr), "Ensembl IDs;", sum(missing_symbol), "had no HGNC symbol;",
    "duplicate symbols collapsed to highest-mean-expression row)\n")

# Sanity check: the 13 DEGs must be present in the normalized matrix
degs <- readRDS("Results_rnaseq/dif_exp_genesymbol_sig.rds")
stopifnot(all(rownames(degs) %in% rownames(normExpr_genesymbol)))

end_time <- Sys.time()
cat("Time taken for analysis: ", format(end_time - start_time), "\n")