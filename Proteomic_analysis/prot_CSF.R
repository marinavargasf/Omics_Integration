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

#Dataframe con info sobre samples
#cardiometa_1 <- read.csv(file =  "Targeted olink/proteomics-CSF-PPEA-D03/protein-expression/releases_2023_v4release_1027_proteomics-CSF-PPEA-D03_olink-explore_protein-expression_CSF-PPEA-D03_cardiometabolic.csv")

#Merged dataframe with all data
proteomic_raw <- rbind(cardiometa, inflammatory, neurology, oncology)
# Convert the first column to row names
rownames(proteomic_raw) <- proteomic_raw[[1]]

# Remove the first column from the dataframe
proteomic_raw <- proteomic_raw[, -1]

#To check that rownames are not repeated 
table(isUnique(rownames(proteomic_raw)))

################################################################################
########################### Create metadata object #############################

#View(samples) #es metadata object 
metadata <- readRDS(file = "~/Documents/Parkinson/proteomic/metadata.rds")

samples <- merge(samples_full, metadata[,c(1,11,13,14)]) 

colnames(samples)[colnames(samples)=="case_control_other_at_baseline"] <- "condition" 

samples <- samples %>% dplyr::select(participant_id, sample_id, condition, sex, age_at_baseline, visit_month)

samples$condition <- factor(samples$condition)

rownames(samples) <- samples$sample_id

# To check how many cases and how many controls
samples_filtrado <- samples %>%
  group_by(participant_id) %>%
  slice(1) %>%
  ungroup()


################################################################################
############################## Merged full data ################################

#check normalization!

full_normalized_counts <- as.matrix(proteomic_raw)

#boxplot(proteomic_raw, las=2)
results <- apply(full_normalized_counts, 2, range)  # Ver rango por columna (muestras)

results_1 <- summary(as.vector(full_normalized_counts)) #log2 normalization assumtion



#boxplot(full_normalized_counts, outline = TRUE, las = 2, main = "Distribución de intensidades proteómicas")

#Heatmap of protein expression
#library("ComplexHeatmap")
#Change "-" for "." for samples to match
samples$sample_id <- gsub("-", ".", samples$sample_id)
#
samples <- samples[order(samples[,"condition"], decreasing = T), , drop = FALSE] # order samples
#sample_ordered <- samples$sample_id
#
#proteomic_to_plot <- proteomic_raw[,sample_ordered]
#
#stopifnot(all(samples$sample_id == colnames(proteomic_to_plot)))# Check
#
#colours <-list("Condition"=c('Case' = "#FFD700", 'Control' = "#008000", 'Other' = "#0077CC"))
#colAnn <- HeatmapAnnotation(df = samples[,"condition", drop=F],
#                            which = 'col',
#                            col = colours,
#                            annotation_width = unit(c(1, 4), 'cm'),
#                            gap = unit(1, 'mm'))
#
#library(circlize)
#col_fun = colorRamp2(c(-max(proteomic_to_plot), 0, max(proteomic_to_plot)), c("indianred", "white", "darkblue"))
#
#proteomic_to_plot <- as.matrix(proteomic_to_plot)
#
#
#hmap <- Heatmap(
#  proteomic_to_plot,
#  #column_split = df_complete.ordered$group,
#  cluster_column_slices = TRUE,
#  name = "Protein expression",
#  #col = greenred(75), 
#  show_row_names = T,
#  show_column_names = T,
#  cluster_rows = TRUE,
#  cluster_columns = TRUE,
#  show_column_dend = TRUE,
#  show_row_dend = TRUE,
#  row_dend_reorder = TRUE,
#  #column_dend_reorder = TRUE,
#  #clustering_method_rows = "ward.D2",
#  #clustering_method_columns = "ward.D2",
#  #width = unit(100, "mm"),
#  top_annotation=colAnn,
#  #top_annotation_height=unit(1.0,"cm"), 
#  row_names_gp = gpar(fontsize = 6),
#  heatmap_legend_param = list(title = "Sample Expression activity")
#)
#
#pheatmap(proteomic_to_plot, 
#         cluster_rows = TRUE, 
#         cluster_cols = TRUE, 
#         show_rownames = FALSE,  # Hide row names for clarity
#         show_colnames = FALSE,  # Hide column names if too many
#         color = col_fun, 
#         main = "Heatmap of Proteomic Data")

hist(full_normalized_counts) #already log normalized


design <- model.matrix(~0 + condition + visit_month + sex + age_at_baseline, samples) # 0 elimina intercept

is.fullrank(design)

# Estimar la correlación intra-individuo
#corfit <- duplicateCorrelation(full_normalized_counts, design, block = samples$participant_id)

# Ver correlación estimada
#corfit$consensus.correlation


# Ajuste del modelo con correlación intra-individuo
#fit <- lmFit(full_normalized_counts, design, block = samples$participant_id, correlation = corfit$consensus.correlation)
fit <- lmFit(full_normalized_counts, design)

# Define contrasts Case VS control
contrast_matrix <- makeContrasts(
  Case_vs_Control = conditionCase - conditionControl,
  levels = design
)

# Fit the contrasts
fit1 <- contrasts.fit(fit, contrast_matrix)

# Compute statistics for differential expression
fit1 <- eBayes(fit1)


# Print the top differentially expressed proteins
top_prots_1 <- topTable(fit1, adjust = "fdr", number = nrow(fit1)) 

top_prots_1<- top_prots_1 %>% filter(abs(logFC)>=0.2, P.Value <=0.05)

print(top_prots_1)

saveRDS(top_prots_1, "~/Documents/Omics_Integration/Proteomic_analysis/Results/differential_expression_proteins_CSF.rds")
################################################################################

# top_prots_1 <- readRDS("~/Documents/Omics_Integration/Proteomic_analysis/Results/differential_expression_proteins_CSF.rds")
# Crear una nueva columna 'Category' e inicializarla como 'Unknown'
top_prots_1$Category <- "Unknown"

# Función para asignar categorías
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

# Aplicar la función a la primera columna del dataframe
top_prots_1$Category <- sapply(rownames(top_prots_1), assign_category)

# Verificar el resultado
top_prots_1

results_proteomics <- rownames(top_prots_1)

#seleccionamos solo las proteinas relevantes para el estudio

filtered_proteins <- full_normalized_counts[rownames(full_normalized_counts) %in% results_proteomics, ]

filtered_proteins <- as.data.frame(t(filtered_proteins))

filtered_proteins$sample_id <- rownames(filtered_proteins)

#Hacemos merge de ambos datasets

filtered_proteins_metadata <- merge(filtered_proteins, samples, by.x ="sample_id", by.y="sample_id")

#Remove all external info, only keeping id and visit_month 
filtered_proteins_metadata <- filtered_proteins_metadata %>% dplyr::select(participant_id, visit_month, all_of(results_proteomics))

saveRDS(filtered_proteins_metadata, "~/Documents/Omics_Integration/Proteomic_analysis/Results/filtered_proteins_metadata_CSF.rds")

