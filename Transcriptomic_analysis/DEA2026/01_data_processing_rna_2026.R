################################################################################
# Combined preprocessing: raw counts -> QC -> neutrophil estimation ->
# case/control assembly -> gene filtering -> aligned metadata for DESeq2 design
################################################################################

library(DESeq2)
library(dplyr)
library(tidyr)
library(purrr)
library(edgeR)
library(biomaRt)

setwd("~/Documents/Omics_Integration/Transcriptomic_analysis/DEA2026/")

################################################################################
## 1. Read metadata, QC, and raw counts
################################################################################

demographics <- read.csv("AMP_PD_Data_v3_2022/clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_Demographics.csv")
QCRNA <- read.csv("AMP_PD_Data_v3_2022/transcriptomics_v4/reports/releases_2023_v4release_1027_rnaseq-WB-RWTS_sequencing_metrics_rna_quality_metrics.csv", row.names = 1)
QCRNA$participant_id <- sub("^([^-]+-[^-]+)-.*", "\\1", rownames(QCRNA))

QCTable <- read.table("AMP_PD_Data_v3_2022/transcriptomics_v4/reports/releases_2023_v4release_1027_rnaseq-WB-RWTS_reports_multiqc_multiqc_data_multiqc_star.txt",
                      sep = "\t", header = TRUE)
QCTable$participant_id <- sub("^([^-]+-[^-]+)-.*", "\\1", QCTable$Sample)
rownames(QCTable) <- QCTable$Sample

demographics <- demographics[demographics$participant_id %in% QCRNA$participant_id, ]
QCRNA <- QCRNA[, grep("Ratio", colnames(QCRNA), invert = TRUE)]
QCTable <- QCTable[QCTable$unmapped_mismatches_percent < 0.05, ]

covariates <- QCRNA %>%
  left_join(demographics %>% dplyr::select(participant_id, sex, age_at_baseline, ethnicity, race),
            by = "participant_id")
rownames(covariates) <- rownames(QCRNA)
covariates <- covariates[rownames(QCTable), ]
covariates <- covariates[covariates$RIN_Value >= 6.5, ]

# guard: left_join can silently create duplicate rows if participant_id isn't 1:1
stopifnot(nrow(covariates) == length(unique(rownames(covariates))))

rawCounts <- read.table("AMP_PD_Data_v3_2022/transcriptomics_v4/subread_feature-counts/releases_2023_v4release_1027_rnaseq-WB-RWTS_subread_feature-counts_matrix.featureCounts.tsv",
                        row.names = 1)
colnames(rawCounts) <- rawCounts[1, ]
rawCounts <- rawCounts[-1, ]
rawCounts <- rawCounts[, colnames(rawCounts) %in% rownames(covariates)]
geneNames <- rownames(rawCounts)
rawCounts <- apply(rawCounts, 2, as.numeric)
rownames(rawCounts) <- sub("\\..*", "", geneNames)
rawCounts <- as.matrix(rawCounts)
covariates <- covariates[colnames(rawCounts), ]

stopifnot(all(colnames(rawCounts) == rownames(covariates)))

################################################################################
## 2. Neutrophil estimation (from hematology labs + 31-gene signature)
################################################################################

hemato <- read.csv("Blood_Chemistry___Hematology_10Mar2025.csv", stringsAsFactors = FALSE)

hemato <- hemato %>%
  filter(EVENT_ID %in% c("BL", "V02", "V04", "V05", "V06", "V08")) %>%
  filter(LTSTNAME != "Which visit being performed?") %>%
  filter(LSIRES != "") %>%
  filter(!LTSTCODE %in% c("HMT70", "HMT71", "HMT98")) %>%
  filter(!is.na(LSIRES)) %>%
  mutate(
    participant_id = paste0("PP-", PATNO),
    Sample = case_when(
      EVENT_ID == "BL"  ~ paste0("PP-", PATNO, "-BLM0T1"),
      EVENT_ID == "V02" ~ paste0("PP-", PATNO, "-SVM6T1"),
      EVENT_ID == "V04" ~ paste0("PP-", PATNO, "-SVM12T1"),
      EVENT_ID == "V05" ~ paste0("PP-", PATNO, "-SVM18T1"),
      EVENT_ID == "V06" ~ paste0("PP-", PATNO, "-SVM24T1"),
      EVENT_ID == "V08" ~ paste0("PP-", PATNO, "-SVM36T1"),
      TRUE ~ NA_character_
    )
  )
hemato$LSIRES <- as.numeric(hemato$LSIRES)

tables <- split(hemato, interaction(hemato$LTSTNAME, hemato$LSIUNIT, drop = TRUE))
cellTypes <- data.frame(Sample = colnames(rawCounts))

for (tablaN in seq_along(tables)) {
  rownames(tables[[tablaN]]) <- tables[[tablaN]]$Sample
  smallD <- data.frame(tables[[tablaN]]$Sample, tables[[tablaN]]$LSIRES)
  colnames(smallD) <- c("Sample", names(table(tables[[tablaN]]$LTSTNAME))[1])
  rownames(smallD) <- smallD$Sample
  cellTypes <- merge(cellTypes, smallD, by = "Sample", all.x = TRUE)
}

cellTypes <- cellTypes[startsWith(cellTypes$Sample, "PP-"), ]
rownames(cellTypes) <- cellTypes$Sample
cellTypes <- cellTypes[, -1]
cellTypes <- cellTypes %>% filter(rowSums(is.na(.)) < ncol(.))

metadataHemato <- merge(covariates, cellTypes, by = "row.names")
rownames(metadataHemato) <- metadataHemato$Row.names
metadataHemato <- metadataHemato[, -1]

metadataHemato <- metadataHemato %>%
  mutate(across(where(~ all(suppressWarnings(!is.na(as.numeric(.))))), as.numeric))
for (variable in colnames(metadataHemato)) {
  if (class(metadataHemato[, variable]) == "character") {
    metadataHemato[, variable] <- factor(metadataHemato[, variable])
  }
}

lmgenes <- c("ENSG00000205927","ENSG00000113389","ENSG00000164047","ENSG00000167186",
             "ENSG00000151012","ENSG00000179299","ENSG00000143226","ENSG00000133687",
             "ENSG00000134827","ENSG00000172292","ENSG00000038219","ENSG00000158488",
             "ENSG00000135373","ENSG00000204345","ENSG00000196549","ENSG00000062716",
             "ENSG00000050748","ENSG00000164187","ENSG00000124721","ENSG00000108846",
             "ENSG00000150051","ENSG00000149516","ENSG00000163162","ENSG00000175793",
             "ENSG00000083290","ENSG00000092969","ENSG00000148572","ENSG00000197208",
             "ENSG00000128594","ENSG00000070731","ENSG00000165046")

ddsVst <- DESeqDataSetFromMatrix(countData = rawCounts, colData = covariates, design = ~1)
ddsVst <- estimateSizeFactors(ddsVst)
vsd <- vst(ddsVst, blind = TRUE)
vstDev <- as.data.frame(t(assay(vsd)))

vstSub <- vstDev[rownames(metadataHemato), lmgenes]
vstSub2 <- merge(vstSub, metadataHemato, by = "row.names")
rownames(vstSub2) <- vstSub2$Row.names
vstSub2 <- vstSub2[, -1]
vstSub2$neutPer <- vstSub2$`Neutrophils (%)`
vstSub2 <- vstSub2[, c(colnames(vstSub), "neutPer")]
vstSub2 <- na.omit(as.data.frame(vstSub2))

fit <- lm(neutPer ~ ., data = vstSub2)
cat("Neutrophil imputation model adj. R^2:", summary(fit)$adj.r.squared, "\n")

neutPredicted <- predict(fit, vstDev)
cat("Predicted range:", range(neutPredicted), "| any out of [0,100]:",
    any(neutPredicted < 0 | neutPredicted > 100), "\n")

covariates$Neutrophils <- neutPredicted
covariates[rownames(vstSub2), "Neutrophils"] <- vstSub2$neutPer

cat("Real neutrophil values used:", sum(rownames(covariates) %in% rownames(vstSub2)),
    "| Imputed:", sum(!rownames(covariates) %in% rownames(vstSub2)), "\n")

stopifnot(!anyNA(covariates$Neutrophils))
attr(covariates, "generated_on") <- Sys.time()
attr(covariates, "script") <- "01_preprocessing_combined.R"
dir.create("./data", recursive = TRUE, showWarnings = FALSE)
saveRDS(covariates, file = "./data/covariates_full.rds")

################################################################################
## 3. Case / control definition
################################################################################

PDHistory <- read.csv("AMP_PD_Data_v3_2022/clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_PD_Medical_History.csv")
PDHistory <- PDHistory[PDHistory$diagnosis != "", ]
number_diagnoses <- rowSums(table(PDHistory$participant_id, PDHistory$diagnosis) > 0)
patientsUniqueDiagnoses <- names(number_diagnoses)[number_diagnoses == 1]
length(patientsUniqueDiagnoses)

clinic <- read.csv("AMP_PD_Data_v3_2022/clinical_metadata_v4_2023/releases_2023_v4release_1027_amp_pd_case_control.csv", row.names = 1)
others <- rownames(clinic)[which(clinic$case_control_other_at_baseline == "Other" | clinic$case_control_other_latest == "Other")]
clinic <- clinic[!rownames(clinic) %in% others, ]
clinic <- clinic[rownames(clinic) %in% patientsUniqueDiagnoses, ]

patients <- clinic[which(clinic$diagnosis_at_baseline %in% c("Idiopathic PD", "Parkinson's Disease")), ]
controls <- clinic[which(clinic$diagnosis_at_baseline == "No PD Nor Other Neurological Disorder"), ]

match_ids <- function(id_list, cols) {
  out <- c()
  for (ptid in id_list) {
    m <- grep(paste0(ptid, "-"), cols, value = TRUE)
    if (length(m) != 0) out <- append(out, m)
  }
  out
}
cases_list    <- match_ids(rownames(patients), colnames(rawCounts))
controls_list <- match_ids(rownames(controls), colnames(rawCounts))

rawCounts <- rawCounts[, colnames(rawCounts) %in% c(cases_list, controls_list)]
dim(rawCounts)

################################################################################
## 4. Gene filtering (low expression + MT/ribo/globin, pinned Ensembl v100)
################################################################################

dge <- DGEList(counts = rawCounts)
keep <- rowSums(cpm(dge) > 10) >= 0.8 * ncol(dge)
dge <- dge[keep, keep.lib.sizes = FALSE]

# raw filtered counts checkpoint (pre-TMM, pre-logCPM) -- this is the DESeq2 input
rawCounts_filtered <- dge$counts

# Pinned to v100 to match GENCODE v29 annotation used for quantification --
# do NOT leave this unpinned: unpinned = today's live release = drifting gene_map
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl") #, version = 100)
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(rawCounts_filtered),
  mart = ensembl
)
mt_genes <- grep("^MT-", gene_map$hgnc_symbol, value = TRUE)
ribo_genes <- grep("^RPL|^RPS", gene_map$hgnc_symbol, value = TRUE)
globin_genes <- grep("HBA|HBB|HBD|HBG", gene_map$hgnc_symbol, value = TRUE)
exclude_genes <- unique(c(ribo_genes, mt_genes, globin_genes))
exclude_genes_ENSEMBL <- gene_map[gene_map$hgnc_symbol %in% exclude_genes, 1]
saveRDS(exclude_genes_ENSEMBL, file = "./data/exclude_genes_ENSEMBL_v100.rds")

rawCounts_filtered <- rawCounts_filtered[!rownames(rawCounts_filtered) %in% exclude_genes_ENSEMBL, ]

################################################################################
## 5. Assemble aligned metadata (covariates + case_control + Cohort + rare-combo filter)
################################################################################

confoundingDataFrame <- covariates[colnames(rawCounts_filtered), ]
confoundingDataFrame <- confoundingDataFrame %>%
  mutate(Cohort = sub("-.*", "", rownames(.)))
confoundingDataFrame$Cohort <- factor(confoundingDataFrame$Cohort)
levels(confoundingDataFrame$Cohort) <- c("BioFIND", "PDBP", "PPMI")
confoundingDataFrame$Position <- factor(confoundingDataFrame$Position)
confoundingDataFrame <- confoundingDataFrame[, !(colnames(confoundingDataFrame) %in% c("participant_id", "visit_name"))]

for (variable in colnames(confoundingDataFrame)) {
  if (class(confoundingDataFrame[, variable]) == "character") {
    confoundingDataFrame[, variable] <- factor(confoundingDataFrame[, variable])
  }
}

confoundingDataFrame$case_control <- NA
confoundingDataFrame[rownames(confoundingDataFrame) %in% cases_list, "case_control"] <- "case"
confoundingDataFrame[rownames(confoundingDataFrame) %in% controls_list, "case_control"] <- "control"
confoundingDataFrame$case_control <- factor(confoundingDataFrame$case_control)

# remove rare Plate x ethnicity combinations (can't be modeled reliably)
combo_freq <- table(confoundingDataFrame$Plate, confoundingDataFrame$ethnicity)
rare_combos <- which(combo_freq == 1, arr.ind = TRUE)
rare_df <- data.frame(
  plate = rownames(combo_freq)[rare_combos[, "row"]],
  ethnicity = colnames(combo_freq)[rare_combos[, "col"]],
  stringsAsFactors = FALSE
)
rows_to_remove <- apply(confoundingDataFrame[, c("Plate", "ethnicity")], 1, function(x) {
  any(x[1] == rare_df$plate & x[2] == rare_df$ethnicity)
})
confoundingDataFrame2 <- confoundingDataFrame[!rows_to_remove, ]

# add visit_month for optional use in the design (pooled-timepoint models)
confoundingDataFrame2$time <- sub("^.*?-.*?-", "", rownames(confoundingDataFrame2))
confoundingDataFrame2$visit_month <- recode(confoundingDataFrame2$time,
                                            "BLM0T1" = 0, "BLM0T1.1" = 0, "SVM0_5T1" = 0,
                                            "SVM6T1" = 6, "SVM6T1.1" = 6,
                                            "SVM12T1" = 12, "SVM12T1.1" = 12,
                                            "SVM18T1" = 18, "SVM18T1.1" = 18,
                                            "SVM24T1" = 24, "SVM24T1.1" = 24,
                                            "SVM36T1" = 36, "SVM36T1.1" = 36
)
unmatched_time <- unique(confoundingDataFrame2$time[is.na(confoundingDataFrame2$visit_month)])
if (length(unmatched_time) > 0) {
  warning("Unmatched 'time' codes found, these samples have NA visit_month: ",
          paste(unmatched_time, collapse = ", "))
}

################################################################################
## 6. Align counts <-> metadata, and save DESeq2-ready checkpoints
################################################################################

common <- intersect(colnames(rawCounts_filtered), rownames(confoundingDataFrame2))
stopifnot(length(common) > 0)
rawCounts_filtered <- rawCounts_filtered[, common]
confoundingDataFrame2 <- confoundingDataFrame2[common, ]
stopifnot(all(colnames(rawCounts_filtered) == rownames(confoundingDataFrame2)))

# make sure no NAs sneak into columns you plan to put in the design formula --
# model.matrix() silently drops NA rows otherwise
design_vars <- c("Neutrophils", "RIN_Value", "Plate", "age_at_baseline",
                 "sex", "ethnicity", "case_control")
na_check <- sapply(confoundingDataFrame2[, design_vars], function(x) sum(is.na(x)))
print(na_check)
stopifnot(sum(na_check) == 0)

saveRDS(rawCounts_filtered, file = "./data/rawCounts_filtered.rds")
saveRDS(confoundingDataFrame2, file = "./data/confoundingDataFrame2_full.rds")

cat("Done. Samples:", ncol(rawCounts_filtered), "| Genes:", nrow(rawCounts_filtered), "\n")
