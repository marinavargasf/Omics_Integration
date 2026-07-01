#########################################################################################################
#Script para la creación de la matriz de counts normalizada y archivo de metadata con sample info y time#
#########################################################################################################

start_time <- Sys.time()

cat("Start time is: ", start_time,"\n")

library("limma")
library("DESeq2")
library("dplyr")
library("edgeR")
library("biomaRt")

#metadata_clinical <- readRDS(file = "metadata_clinical_variables.rds")

countData <- readRDS(file = "expression_crude.rds")
case_control <- readRDS(file = "casecontrolDF.rds")
covariates <- readRDS(file = "confoundingDF.rds")
neutEst <- readRDS(file = "neutPer_casecontrol.rds")

rownames(countData) <- sub("\\..*", "", rownames(countData))


keep <- rowSums(cpm(countData) > 10) >= 0.8 * ncol(countData)
countData <- countData[keep,]

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

gene_map <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = rownames(countData),
    mart = ensembl
)

mt_genes <- grep("^MT-", gene_map$hgnc_symbol, value=TRUE)
ribo_genes <- grep("^RPL|^RPS", gene_map$hgnc_symbol, value=TRUE)
globin_genes <- grep("HBA|HBB|HBD|HBG", gene_map$hgnc_symbol, value=TRUE)
exclude_genes <- unique(c(ribo_genes, mt_genes, globin_genes))
exclude_genes_ENSEMBL <- gene_map[gene_map$hgnc_symbol %in% exclude_genes, 1]

countData <-countData[!rownames(countData) %in% exclude_genes_ENSEMBL,]

#################### metadata object #############################

metadata <- case_control

rownames(metadata) <- metadata$sample_id

metadata <- merge(metadata, neutEst, by="row.names", all = T)

rownames(metadata) <- metadata$Row.names
metadata$Row.names <- NULL

colnames(metadata)[colnames(metadata) == "status"] <- "condition"

metadata$participant_id <- sub("^((.*?)-(.*?))-.*", "\\1", rownames(metadata))
metadata$time <- sub("^.*?-.*?-", "", rownames(metadata))

metadata$visit_month <- recode(metadata$time,
                                       "BLM0T1" = 0,
                                       "BLM0T1.1" = 0,
                                       "SVM0_5T1" = 0,
                                       "SVM12T1" = 12,
                                       "SVM12T1.1" = 12,
                                       "SVM18T1" = 18,
                                       "SVM24T1" = 24,
                                       "SVM24T1.1" = 24,
                                       "SVM36T1" = 36,
                                       "SVM36T1.1" = 36,
                                       "SVM6T1" = 6,
                                       "SVM6T1.1" = 6,
)

metadata <- merge(metadata, covariates, by="row.names") #nos quedamos solo con las muestras en metadata (ya filtrados)

rownames(metadata) <- metadata$Row.names
metadata$Row.names <- NULL
metadata$condition <- factor(metadata$condition)


#Filtramos countData para que tenga las mismas samples que metadata! error feb2025 un sample no tenia covariates!! design con NA
common_samples <- intersect(colnames(countData), rownames(metadata))
countData <- countData[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

#Normalize data###################################

dge <- DGEList(counts= countData)

normFactors = calcNormFactors(countData, method="TMM", refColumn = 1, logratioTrim = 0.3, sumTrim = 0.05, doWeighting = T, Acutoff = -1e+10)
normFactors = normFactors*(colSums(countData)/mean(colSums(countData)))
names(normFactors) = colnames(countData)
#countData <- round(countData,0)

logCPM <- cpm(dge, log = T , prior.count = 1)


selectedVars=c("neutPer","RIN_Value","Submitted_Volume__ul_","Plate","Concentration","sex","age_at_baseline","race")

matrix2Correct <- metadata[, selectedVars]
design <- model.matrix(~., data=matrix2Correct)[,-1]

library(limma)
normExpr <- removeBatchEffect(logCPM,
                              covariates = design)


saveRDS(metadata, file = "Normalization/metadata.rds")
saveRDS(normExpr, file = "Normalization/normalized_expression.rds")
#save.image("processed_data.Data")


cat("Normalization done \n")

#print(is.fullrank(design))


#Print time
end_time <- Sys.time()

time_taken <- end_time - start_time

cat("Time taken for analysis: ", time_taken, "n")
