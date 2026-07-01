################################################################################
################################################################################

#Record the start time
start_time <- Sys.time()

cat("Start time is: ", start_time,"\n")

library("DESeq2")
library("dplyr")
library("edgeR")
library("biomaRt")

################################################################################
################################## Read data ###################################

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

################################################################################
########################### Create metadata object #############################

#eliminamos NAS

#covariates <- na.omit(covariates) 

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
#metadata$condition <- factor(metadata$condition)
#saveRDS(metadata, file = "Results_feb25/metadata.rds")

head(metadata)

cat("Original metadata \n")


#Check to see differences in samples

#Filtramos countData para que tenga las mismas samples que metadata! error feb2025 un sample no tenia covariates!! design con NA
common_samples <- intersect(colnames(countData), rownames(metadata))
countData <- countData[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]


# Comparar si los nombres son exactamente iguales
#Reordenamos countData

countData <-countData[, rownames(metadata)] 
cat("ya reordenado!")
#head(countData)
################################################################################
#############Differential expression on rnaseq data#############################

normFactors = calcNormFactors(countData, method="TMM", refColumn = 1, logratioTrim = 0.3, sumTrim = 0.05, doWeighting = T, Acutoff = -1e+10)
normFactors = normFactors*(colSums(countData)/mean(colSums(countData)))
names(normFactors) = colnames(countData)
#countData <- round(countData,0)

cat("Normalization factors are calculated \n")

design <- model.matrix(~condition + visit_month + sex + age_at_baseline + race + neutPer + Plate + RIN_Value + Submitted_Volume__ul_  + X_260_280_Ratio, metadata)

cat("Is design fullrank? \n")
print(is.fullrank(design))

cat("Does metadata and countdata match? must be TRUE")
print(all(rownames(metadata)==colnames(countData)))

cat("Metadata dim \n")
print(dim(metadata))

cat("countData dim \n")
print(dim(countData))

cat("Ultimo check \n")
print(ncol(design))

ddsFromMatrix <- DESeqDataSetFromMatrix(countData = countData,
                                        colData = metadata,
                                        design = design) 

 
# Desactiva la normalización por factores de tamaño, para cuando tenemos datos ya normalizados 
#sizeFactors(ddsFromMatrix) <- rep(1, ncol(ddsFromMatrix)) 

cat("Ahora empieza DESeq\n")

sizeFactors(ddsFromMatrix) = normFactors
dds <- DESeq(ddsFromMatrix, parallel = T)

cat("DESeq done!\n")

save.image("Results_junio/deseq2_output.RData")


#res <- as.data.frame(results(dds, name = c("conditionControl")))

#save.image("Results_v3/test_final_long.RData")

#saveRDS(res, file = "Results_mayo25/res_deseq.rds")


#Print time
end_time <- Sys.time()

time_taken <- end_time - start_time

cat("Time taken for analysis: ", time_taken, "\n")
