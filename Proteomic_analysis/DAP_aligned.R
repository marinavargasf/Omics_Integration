################################################################################
# CSF Proteomics — Differential Expression Analysis (Case vs Control)
################################################################################

setwd("~/Documents/Omics_Integration/Proteomic_analysis/")

library("limma")
library("dplyr")
library("edgeR")
library("Biobase")

samples_full <- read.csv(file="~/Documents/Parkinson/proteomic/Targeted olink/proteomics-CSF-PPEA-D03/releases_2023_v4release_1027_proteomics-CSF-PPEA-D03_proteomics_CSF-PPEA-D03_samples.csv")
glossary <- read.csv(file="~/Documents/Parkinson/proteomic/Targeted olink/proteomics-CSF-PPEA-D03/releases_2023_v4release_1027_proteomics-CSF-PPEA-D03_ColumnHeaderGlossary.csv")

cardiometa <- read.csv(file =  "~/Documents/Parkinson/proteomic/Targeted olink/proteomics-CSF-PPEA-D03/protein-expression/releases_2023_v4release_1027_proteomics-CSF-PPEA-D03_olink-explore_protein-expression_CSF-PPEA-D03_matrix_cardiometabolic.csv")
inflammatory <- read.csv(file =  "~/Documents/Parkinson/proteomic/Targeted olink/proteomics-CSF-PPEA-D03/protein-expression/releases_2023_v4release_1027_proteomics-CSF-PPEA-D03_olink-explore_protein-expression_CSF-PPEA-D03_matrix_inflammation.csv")
neurology <- read.csv(file =  "~/Documents/Parkinson/proteomic/Targeted olink/proteomics-CSF-PPEA-D03/protein-expression/releases_2023_v4release_1027_proteomics-CSF-PPEA-D03_olink-explore_protein-expression_CSF-PPEA-D03_matrix_neurology.csv")
oncology <- read.csv(file =  "~/Documents/Parkinson/proteomic/Targeted olink/proteomics-CSF-PPEA-D03/protein-expression/releases_2023_v4release_1027_proteomics-CSF-PPEA-D03_olink-explore_protein-expression_CSF-PPEA-D03_matrix_oncology.csv")

#Merged dataframe with all data
proteomic_raw <- rbind(cardiometa, inflammatory, neurology, oncology)
rownames(proteomic_raw) <- proteomic_raw[[1]]
proteomic_raw <- proteomic_raw[, -1]

table(isUnique(rownames(proteomic_raw)))

################################################################################
########################### Create metadata object #############################

metadata <- readRDS(file = "~/Documents/Parkinson/proteomic/metadata.rds")

samples <- merge(samples_full, metadata[,c(1,11,13,14)]) 

colnames(samples)[colnames(samples)=="case_control_other_at_baseline"] <- "condition" 

samples <- samples %>% dplyr::select(participant_id, sample_id, condition, sex, age_at_baseline, visit_month)

samples$condition <- factor(samples$condition)

rownames(samples) <- samples$sample_id

samples_filtrado <- samples %>%
  group_by(participant_id) %>%
  dplyr::slice(1) %>%
  ungroup()

################################################################################
############################## Merged full data ################################

full_normalized_counts <- as.matrix(proteomic_raw)

results <- apply(full_normalized_counts, 2, range)
results_1 <- summary(as.vector(full_normalized_counts)) #log2 normalization assumtion

#Change "-" for "." for samples to match
samples$sample_id <- gsub("-", ".", samples$sample_id)
#
samples <- samples[order(samples[,"condition"], decreasing = T), , drop = FALSE] # order samples
table(samples$condition)
# Exclude "Other" diagnosis before fitting -- residual variance and shared
# covariate coefficients (visit_month, sex, age) should not be estimated
# using a group with unverified covariate structure
samples <- samples %>% filter(condition %in% c("Case", "Control"))
samples$condition <- droplevels(samples$condition)

# Full sample set (Case/Control only), explicitly aligned by name
full_normalized_counts_full <- full_normalized_counts[, samples$sample_id]
stopifnot(all(samples$sample_id == colnames(full_normalized_counts_full)))

hist(full_normalized_counts_full) #already log normalized

design <- model.matrix(~0 + condition + visit_month + sex + age_at_baseline, samples) # 0 elimina intercept

is.fullrank(design)

fit <- lmFit(full_normalized_counts_full, design)

# Define contrasts Case VS control
contrast_matrix <- makeContrasts(
  Case_vs_Control = conditionCase - conditionControl,
  levels = design
)

fit1 <- contrasts.fit(fit, contrast_matrix)
fit1 <- eBayes(fit1)

top_prots_full <- topTable(fit1, adjust = "fdr", number = nrow(fit1))
top_prots_sig <- top_prots_full %>% filter(abs(logFC)>=0.2, adj.P.Val <=0.05)

print(top_prots_sig)
saveRDS(top_prots_full, "~/Documents/Omics_Integration/Proteomic_analysis/Results_full_aligned/full_differential_expression_proteins_CSF.rds")
saveRDS(top_prots_sig, "~/Documents/Omics_Integration/Proteomic_analysis/Results_full_aligned/differential_expression_proteins_CSF.rds")
write.csv(top_prots_sig, "~/Documents/Omics_Integration/Proteomic_analysis/Results_full_aligned/differential_expression_proteins_CSF.csv", row.names = TRUE)

################################################################################

top_prots_sig$Category <- "Unknown"

assign_category <- function(protein) {
  if (protein %in% cardiometa$UniProt) {
    return("Cardiometa")
  } else if (protein %in% inflammatory$UniProt) {
    return("Inflammatory")
  } else if (protein %in% neurology$UniProt) {
    return("Neurology")
  } else if (protein %in% oncology$UniProt) {
    return("Oncology")
  } else {
    return("Unknown")
  }
}

top_prots_sig$Category <- sapply(rownames(top_prots_sig), assign_category)

top_prots_sig

results_proteomics <- rownames(top_prots_sig)

filtered_proteins <- full_normalized_counts_full[rownames(full_normalized_counts_full) %in% results_proteomics, ]

filtered_proteins <- as.data.frame(t(filtered_proteins))

filtered_proteins$sample_id <- rownames(filtered_proteins)

filtered_proteins_metadata <- merge(filtered_proteins, samples, by.x = "sample_id", by.y = "sample_id")

filtered_proteins_metadata <- filtered_proteins_metadata %>% dplyr::select(participant_id, visit_month, all_of(results_proteomics))

saveRDS(filtered_proteins_metadata, "~/Documents/Omics_Integration/Proteomic_analysis/Results_full_aligned/filtered_proteins_metadata_CSF.rds")



################################################################################
# Plasma Proteomics — Differential Expression Analysis (Case vs Control)
################################################################################

setwd("~/Documents/Omics_Integration/Proteomic_analysis/")

library("limma")
library("dplyr")
library("edgeR")
library("Biobase")

samples_full <- read.csv(file="~/Documents/Parkinson/proteomic/Targeted olink/proteomics-PLA-PPEA-D03/releases_2023_v4release_1027_proteomics-PLA-PPEA-D03_proteomics_PLA-PPEA-D03_samples.csv")
glossary <- read.csv(file="~/Documents/Parkinson/proteomic/Targeted olink/proteomics-PLA-PPEA-D03/releases_2023_v4release_1027_proteomics-PLA-PPEA-D03_ColumnHeaderGlossary.csv")

cardiometa <- read.csv(file =  "~/Documents/Parkinson/proteomic/Targeted olink/proteomics-PLA-PPEA-D03/protein-expression/releases_2023_v4release_1027_proteomics-PLA-PPEA-D03_olink-explore_protein-expression_PLA-PPEA-D03_matrix_cardiometabolic.csv")
inflammatory <- read.csv(file =  "~/Documents/Parkinson/proteomic/Targeted olink/proteomics-PLA-PPEA-D03/protein-expression/releases_2023_v4release_1027_proteomics-PLA-PPEA-D03_olink-explore_protein-expression_PLA-PPEA-D03_matrix_inflammation.csv")
neurology <- read.csv(file =  "~/Documents/Parkinson/proteomic/Targeted olink/proteomics-PLA-PPEA-D03/protein-expression/releases_2023_v4release_1027_proteomics-PLA-PPEA-D03_olink-explore_protein-expression_PLA-PPEA-D03_matrix_neurology.csv")
oncology <- read.csv(file =  "~/Documents/Parkinson/proteomic/Targeted olink/proteomics-PLA-PPEA-D03/protein-expression/releases_2023_v4release_1027_proteomics-PLA-PPEA-D03_olink-explore_protein-expression_PLA-PPEA-D03_matrix_oncology.csv")

proteomic_raw <- rbind(cardiometa, inflammatory, neurology, oncology)
rownames(proteomic_raw) <- proteomic_raw[[1]]
proteomic_raw <- proteomic_raw[, -1]

table(isUnique(rownames(proteomic_raw)))

################################################################################
########################### Create metadata object ############################

metadata <- readRDS(file = "~/Documents/Parkinson/proteomic/metadata.rds")

samples <- merge(samples_full, metadata[,c(1,11,13,14)])

colnames(samples)[colnames(samples)=="case_control_other_at_baseline"] <- "condition"

samples <- samples %>% dplyr::select(participant_id, sample_id, condition, sex, age_at_baseline, visit_month)

samples$condition <- factor(samples$condition)

rownames(samples) <- samples$sample_id

samples_filtrado <- samples %>%
  group_by(participant_id) %>%
  dplyr::slice(1) %>%
  ungroup()

################################################################################
############################## Merged full data #################################

full_normalized_counts <- as.matrix(proteomic_raw)

results <- apply(full_normalized_counts, 2, range)
results_1 <- summary(as.vector(full_normalized_counts)) #log2 normalization assumtion

#Change "-" for "." for samples to match
samples$sample_id <- gsub("-", ".", samples$sample_id)
#
samples <- samples[order(samples[,"condition"], decreasing = T), , drop = FALSE] # order samples
table(samples$condition)
# Exclude "Other" diagnosis before fitting (same rationale as CSF above)
samples <- samples %>% filter(condition %in% c("Case", "Control"))
samples$condition <- droplevels(samples$condition)

full_normalized_counts_full <- full_normalized_counts[, samples$sample_id]
stopifnot(all(samples$sample_id == colnames(full_normalized_counts_full)))

hist(full_normalized_counts_full) #already log normalized

design <- model.matrix(~0 + condition + visit_month + sex + age_at_baseline, samples) # 0 elimina intercept

is.fullrank(design)

fit <- lmFit(full_normalized_counts_full, design)

contrast_matrix <- makeContrasts(
  Case_vs_Control = conditionCase - conditionControl,
  levels = design
)

fit1 <- contrasts.fit(fit, contrast_matrix)
fit1 <- eBayes(fit1)

top_prots_full <- topTable(fit1, adjust = "fdr", number = nrow(fit1))
top_prots_sig <- top_prots_full %>% filter(abs(logFC)>=0.2, adj.P.Val <=0.05)

print(top_prots_sig)
saveRDS(top_prots_full, "~/Documents/Omics_Integration/Proteomic_analysis/Results_full_aligned/full_differential_expression_proteins_PLA.rds")
saveRDS(top_prots_sig, "~/Documents/Omics_Integration/Proteomic_analysis/Results_full_aligned/differential_expression_proteins_PLA.rds")
write.csv(top_prots_sig, "~/Documents/Omics_Integration/Proteomic_analysis/Results_full_aligned/differential_expression_proteins_PLA.csv", row.names = TRUE)

################################################################################

top_prots_sig$Category <- "Unknown"

assign_category <- function(protein) {
  if (protein %in% cardiometa$UniProt) {
    return("Cardiometa")
  } else if (protein %in% inflammatory$UniProt) {
    return("Inflammatory")
  } else if (protein %in% neurology$UniProt) {
    return("Neurology")
  } else if (protein %in% oncology$UniProt) {
    return("Oncology")
  } else {
    return("Unknown")
  }
}

top_prots_sig$Category <- sapply(rownames(top_prots_sig), assign_category)

top_prots_sig

results_proteomics <- rownames(top_prots_sig)

filtered_proteins <- full_normalized_counts_full[rownames(full_normalized_counts_full) %in% results_proteomics, ]

filtered_proteins <- as.data.frame(t(filtered_proteins))

filtered_proteins$sample_id <- rownames(filtered_proteins)

filtered_proteins_metadata <- merge(filtered_proteins, samples, by.x = "sample_id", by.y = "sample_id")

filtered_proteins_metadata <- filtered_proteins_metadata %>% dplyr::select(participant_id, visit_month, all_of(results_proteomics))

saveRDS(filtered_proteins_metadata, "~/Documents/Omics_Integration/Proteomic_analysis/Results_full_aligned/filtered_proteins_metadata_PLA.rds")
