################################################################################
# Combined preprocessing: raw counts -> QC -> neutrophil estimation ->
# case/control assembly -> gene filtering -> aligned metadata for DESeq2 design
# + diagnostic swamp/PCA plots before and after batch correction
################################################################################

library(DESeq2)
library(dplyr)
library(tidyr)
library(purrr)
library(edgeR)
library(biomaRt)
library(swamp)
library(limma)
library(ggplot2)

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
attr(covariates, "script") <- "01_data_processing_rna_2026.R"
dir.create("./data", recursive = TRUE, showWarnings = FALSE)
saveRDS(covariates, file = "./data/covariates_full.rds")

################################################################################
## 3. Case / control definition
################################################################################

PDHistory <- read.csv("AMP_PD_Data_v3_2022/clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_PD_Medical_History.csv")
PDHistory <- PDHistory[PDHistory$diagnosis != "", ]
number_diagnoses <- rowSums(table(PDHistory$participant_id, PDHistory$diagnosis) > 0)
patientsUniqueDiagnoses <- names(number_diagnoses)[number_diagnoses == 1]

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

################################################################################
## 4. Gene filtering (low expression + MT/ribo/globin, pinned Ensembl v100)
################################################################################

dge <- DGEList(counts = rawCounts)
keep <- rowSums(cpm(dge) > 10) >= 0.8 * ncol(dge)
dge <- dge[keep, keep.lib.sizes = FALSE]

rawCounts_filtered <- dge$counts

exclude_map_path <- "./data/exclude_genes_ENSEMBL_v100.rds"

if (file.exists(exclude_map_path)) {
  exclude_genes_ENSEMBL <- readRDS(exclude_map_path)
  cat("Loaded cached MT/ribo/globin exclusion list (", length(exclude_genes_ENSEMBL),
      "genes) from", exclude_map_path, "-- skipping Ensembl query\n")
} else {
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")#, version = 100)
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
  saveRDS(exclude_genes_ENSEMBL, file = exclude_map_path)
}

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

combo_freq <- table(confoundingDataFrame$Plate, confoundingDataFrame$race)
rare_combos <- which(combo_freq == 1, arr.ind = TRUE)
rare_df <- data.frame(
  plate = rownames(combo_freq)[rare_combos[, "row"]],
  race = colnames(combo_freq)[rare_combos[, "col"]],
  stringsAsFactors = FALSE
)
rows_to_remove <- apply(confoundingDataFrame[, c("Plate", "race")], 1, function(x) {
  any(x[1] == rare_df$plate & x[2] == rare_df$race)
})
#confoundingDataFrame2 <- confoundingDataFrame[!rows_to_remove, ]

confoundingDataFrame2 <- confoundingDataFrame
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

design_vars <- c("Neutrophils", "RIN_Value", "Plate", "age_at_baseline",
                 "sex", "race", "case_control", "Submitted_Volume__ul_")
na_check <- sapply(confoundingDataFrame2[, design_vars], function(x) sum(is.na(x)))
print(na_check)
stopifnot(sum(na_check) == 0)

saveRDS(rawCounts_filtered, file = "./data/rawCounts_filtered.rds")
saveRDS(confoundingDataFrame2, file = "./data/confoundingDataFrame2_full.rds")

cat("Done. Samples:", ncol(rawCounts_filtered), "| Genes:", nrow(rawCounts_filtered), "\n")

################################################################################
################################################################################
## 7. DIAGNOSTIC PLOTS: swamp + PCA, before and after batch correction
##    (for exploratory/QC purposes only -- does NOT alter rawCounts_filtered,
##    which is corrected via the DESeq2 design formula, not removeBatchEffect)
################################################################################
################################################################################

dir.create("results/Bias_correction", recursive = TRUE, showWarnings = FALSE)

# ---- Build logCPM for diagnostics (TMM-normalized, log scale) --------------
dge_diag <- DGEList(counts = rawCounts_filtered)
dge_diag <- calcNormFactors(dge_diag, method = "TMM")
logCPM <- cpm(dge_diag, log = TRUE, prior.count = 1)

# align diagnostic metadata to logCPM's sample order
diag_meta <- confoundingDataFrame2[colnames(logCPM), ]
stopifnot(all(colnames(logCPM) == rownames(diag_meta)))
diag_meta[] <- lapply(diag_meta, function(x) if (is.character(x)) factor(x) else x)

# --------------------------------------------------------------------------
## 7a. swamp::prince() BEFORE correction -- heatmap of variance explained
##     per principal component by each covariate
# --------------------------------------------------------------------------

resConf_before <- prince(logCPM, diag_meta, top = 10)
saveRDS(resConf_before, "results/Bias_correction/resConf_before.rds")
#resConf_before <- readRDS("results/Bias_correction/resConf_before.rds")

png("results/Bias_correction/swamp_before_correction.png", width = 960, height = 960)
prince.plot(resConf_before, margins = c(5, 15))
dev.off()

pdf("results/Bias_correction/swamp_before_correction.pdf", width = 10, height = 10)
prince.plot(resConf_before, margins = c(5, 15))
dev.off()

# --------------------------------------------------------------------------
## 7b. PCA scatter plots BEFORE correction -- PC pair per variable chosen
##     automatically from resConf_before$linp (the swamp/prince() p-value
##     matrix), instead of hardcoded PC indices.
# --------------------------------------------------------------------------

pca_before <- prcomp(t(logCPM), scale. = TRUE)
pca_before_df <- as.data.frame(pca_before$x[, 1:10])
pca_before_df <- cbind(pca_before_df, diag_meta[rownames(pca_before_df), ])
var_importance_before <- summary(pca_before)$importance[2, ]

# Pick, for each variable (row) in linp, the two PCs (columns) with the
# smallest p-value -- i.e. the PCs where that variable actually showed signal.
# NOTE: resConf_before$linp has NO column names (colnames are NULL) -- columns
# are positional, column i corresponds to PCi. Do NOT try to parse PC numbers
# out of names() here, that silently returns integer(0) for every row.
get_top_pcs <- function(linp_matrix, n_pcs = 2) {
  result <- apply(linp_matrix, 1, function(row) {
    order(row)[1:min(n_pcs, length(row))]   # column positions = PC numbers directly
  }, simplify = FALSE)
  names(result) <- rownames(linp_matrix)   # apply() does not reliably preserve rownames as
                                            # list names -- assign explicitly, or the for-loop
                                            # below silently iterates zero times
  result
}

pc_pairs_before <- get_top_pcs(resConf_before$linp, n_pcs = 2)
cat("pc_pairs_before: ", length(pc_pairs_before), "variables, names:",
    paste(names(pc_pairs_before), collapse = ", "), "\n")

for (v in names(pc_pairs_before)) {
  if (!v %in% colnames(pca_before_df)) next
  pair <- pc_pairs_before[[v]]
  if (length(pair) < 2) next
  pcx <- paste0("PC", pair[1]); pcy <- paste0("PC", pair[2])
  if (!(pcx %in% colnames(pca_before_df)) || !(pcy %in% colnames(pca_before_df))) next

  p <- ggplot(pca_before_df, aes(x = .data[[pcx]], y = .data[[pcy]], color = .data[[v]])) +
    geom_point(alpha = 0.7, size = 1.5) +
    labs(title = paste("PCA before correction —", v, "(", pcx, "vs", pcy, ")"),
         x = paste0(pcx, " (", round(var_importance_before[pair[1]] * 100, 1), "%)"),
         y = paste0(pcy, " (", round(var_importance_before[pair[2]] * 100, 1), "%)")) +
    theme_minimal()

  if (is.numeric(pca_before_df[[v]])) p <- p + scale_color_viridis_c()

  ggsave(paste0("results/Bias_correction/PCA_before_", v, "_", pcx, "vs", pcy, ".png"),
         p, width = 7, height = 6)
}

# --------------------------------------------------------------------------
## 7c. Batch correction with removeBatchEffect (diagnostic version only)
# --------------------------------------------------------------------------

selectedVars <- c("visit_month", "sex", "age_at_baseline", "race",
                   "Neutrophils", "Plate", "RIN_Value", "Submitted_Volume__ul_")
selectedVars <- selectedVars[selectedVars %in% colnames(diag_meta)]

na_check_diag <- sapply(diag_meta[, selectedVars], function(x) sum(is.na(x)))
stopifnot(sum(na_check_diag) == 0)

matrix2Correct <- diag_meta[, selectedVars]
designCovariates <- model.matrix(~., data = matrix2Correct)[, -1]

normExpr_diag <- removeBatchEffect(logCPM, covariates = designCovariates)

# --------------------------------------------------------------------------
## 7d. swamp::prince() AFTER correction
# --------------------------------------------------------------------------

resConf_after <- prince(as.matrix(normExpr_diag), diag_meta[colnames(normExpr_diag), ], top = 10)
saveRDS(resConf_after, "results/Bias_correction/resConf_after.rds")
#resConf_after <- readRDS("results/Bias_correction/resConf_after.rds")

png("results/Bias_correction/swamp_after_correction.png", width = 960, height = 960)
prince.plot(resConf_after, margins = c(5, 15))
dev.off()

pdf("results/Bias_correction/swamp_after_correction.pdf", width = 10, height = 10)
prince.plot(resConf_after, margins = c(5, 15))
dev.off()

# --------------------------------------------------------------------------
## 7e. PCA scatter plots AFTER correction -- PC pair per variable chosen
##     automatically from resConf_after$linp, same logic as 7b. Uses the
##     post-correction significance, so a variable can land on a different
##     PC pair after correction than it did before.
# --------------------------------------------------------------------------

pca_after <- prcomp(t(normExpr_diag), scale. = TRUE)
pca_after_df <- as.data.frame(pca_after$x[, 1:10])
pca_after_df <- cbind(pca_after_df, diag_meta[rownames(pca_after_df), ])
var_importance_after <- summary(pca_after)$importance[2, ]

pc_pairs_after <- get_top_pcs(resConf_after$linp, n_pcs = 2)
cat("pc_pairs_after: ", length(pc_pairs_after), "variables, names:",
    paste(names(pc_pairs_after), collapse = ", "), "\n")

for (v in names(pc_pairs_after)) {
  if (!v %in% colnames(pca_after_df)) next
  pair <- pc_pairs_after[[v]]
  if (length(pair) < 2) next
  pcx <- paste0("PC", pair[1]); pcy <- paste0("PC", pair[2])
  if (!(pcx %in% colnames(pca_after_df)) || !(pcy %in% colnames(pca_after_df))) next

  p <- ggplot(pca_after_df, aes(x = .data[[pcx]], y = .data[[pcy]], color = .data[[v]])) +
    geom_point(alpha = 0.7, size = 1.5) +
    labs(title = paste("PCA after correction —", v, "(", pcx, "vs", pcy, ")"),
         x = paste0(pcx, " (", round(var_importance_after[pair[1]] * 100, 1), "%)"),
         y = paste0(pcy, " (", round(var_importance_after[pair[2]] * 100, 1), "%)")) +
    theme_minimal()

  if (is.numeric(pca_after_df[[v]])) p <- p + scale_color_viridis_c()

  ggsave(paste0("results/Bias_correction/PCA_after_", v, "_", pcx, "vs", pcy, ".png"),
         p, width = 7, height = 6)
}

cat("Diagnostic swamp/PCA plots saved to results/Bias_correction/\n")

resConf_after <- readRDS("~/Documents/Omics_Integration/Transcriptomic_analysis/DEA2026/results/Bias_correction/resConf_after.rds")
resConf_before <- readRDS("~/Documents/Omics_Integration/Transcriptomic_analysis/DEA2026/results/Bias_correction/resConf_before.rds")

pc_labels <- function(res) {
  k <- ncol(res$linp)
  paste0("PC", seq_len(k), " (", round(res$prop[seq_len(k)], 0), "%)")
}

show_prince <- function(res, titulo) {
  p_tab  <- signif(res$linp, 3)
  r2_tab <- round(res$rsquared, 3)
  colnames(p_tab) <- colnames(r2_tab) <- pc_labels(res)
  cat("\n====", titulo, "— p-valores ====\n"); print(p_tab)
  cat("\n====", titulo, "— R² ====\n");       print(r2_tab)
}

show_prince(resConf_before, "Transcriptómica ANTES")
show_prince(resConf_after,  "Transcriptómica DESPUÉS")


#resConfCSF_after <- readRDS("~/Documents/Omics_Integration/Proteomic_analysis/Results_full_aligned/Bias_correction_raw/resConf_proteomics_CSF_after.rds")
#resConfCSF_before <- readRDS("~/Documents/Omics_Integration/Proteomic_analysis/Results_full_aligned/Bias_correction_raw/resConf_proteomics_CSF_before.rds")

#show_prince(resConfCSF_before, "CSF ANTES")
#show_prince(resConfCSF_after,  "CSF DESPUÉS")

#resConfPLA_after <- readRDS("~/Documents/Omics_Integration/Proteomic_analysis/Results_full_aligned/Bias_correction_raw/resConf_proteomics_PLA_after.rds")
#resConfPLA_before <- readRDS("~/Documents/Omics_Integration/Proteomic_analysis/Results_full_aligned/Bias_correction_raw/resConf_proteomics_PLA_before.rds")

#show_prince(resConfPLA_before, "PLA ANTES")
#show_prince(resConfPLA_after,  "PLA DESPUÉS")
