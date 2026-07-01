################################################################################

#Automatized script for longitudinal analysis by linear mixed models

#Marina Vargas-Fernández


library(dplyr)      # filtering datasets
library(rstatix)    # summary statistics
library(ggpubr)     # convenient summary statistics and plots
library(tidyverse)  # not sure if usefull
library(lme4)       # linear mixed effects model
library(qgraph)     # network package
library(MariNET)

setwd("~/Documents/Parkinson/")

#Run required previous scripts

source("clinical_metadata_v4_2023/DataReading_v4.R")
source("shifted_matrix.R")
source("data_processing.R") #previous NetwrokComparison_3.R
#source("Long_analysis_v4.R")

#1. Collect and filter data 
scores <- readRDS("GSEA/pathwayScores_samplesFiltered_GO.rds")
NES <- readRDS("GSEA/nes_GO.RDS")
NES$regulation <- ifelse(
  NES$NES > 0, "Positive regulated",
  ifelse(NES$NES < 0, "Negative regulated", NA)
)
NES$pathway <- gsub("HALLMARK_", "", NES$pathway)


# View the updated dataframe
#head(fgsea_results)

#df_filtered <- data_normalized
#Creation of gene_Expr matrix with sample_id, visit_month
selected_pathways <- as.data.frame(t(scores))
pathway_names <- rownames(scores)

selected_pathways$participant_id <- sub("^((.*?)-(.*?))-.*", "\\1", rownames(selected_pathways))

selected_pathways$time <- sub("^.*?-.*?-", "", rownames(selected_pathways))

selected_pathways$visit_month <- recode(selected_pathways$time,
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

#saveRDS(selected_pathways, "GSEA/selected_pathways_longformat.rds")

#To check how many samples we have
samples_full <- merge(selected_pathways, metadata, by ="participant_id", all.x = T)
# To check how many cases and how many controls
samples_filtrado <- samples_full %>%
  group_by(participant_id) %>%
  slice(1) %>%
  ungroup()

table(samples_filtrado$diagnosis_at_baseline)
# View(metadata)

clinical_data_2 <- merge(clinical_data, selected_pathways, by.x = c("participant_id", "visit_month"), by.y =c("participant_id", "visit_month"), all=TRUE)

master_matrix <- merge(clinical_data_2, metadata,by = "participant_id", all.x = T) 

# left join, only need clinical information on patients who have registers 
dim(master_matrix)

#Filter only cases to check influence on gene expression

df_no_controls <- subset(master_matrix, diagnosis_at_baseline =="Idiopathic PD")


table(df_no_controls$case_control_other_at_baseline)

library(dplyr)

df_no_controls <- df_no_controls %>% dplyr::select( participant_id, visit_month, study, diagnosis_at_baseline, age_at_baseline, sex, race, case_control_other_at_baseline, global_famhistory, has_known_GBA_mutation_in_WGS, has_known_LRRK2_mutation_in_WGS,
                                            has_known_SNCA_mutation_in_WGS, has_known_APOE_E4_mutation_in_WGS, has_known_PD_Mutation_in_WGS, on_levodopa, on_dopamine_agonist, on_other_pd_medications, mds_updrs_part_i_summary_score, mds_updrs_part_ii_summary_score,
                                            mds_updrs_part_iii_summary_score, mds_updrs_part_iv_summary_score, moca_total_score, pdq39_mobility_score, pdq39_adl_score, pdq39_emotional_score, pdq39_stigma_score, pdq39_social_score, pdq39_cognition_score, pdq39_communication_score, pdq39_discomfort_score, 
                                            Schwad_ADL_score, REM_score, ess_sleepiness_score, upsit_total_score,
                                            all_of(pathway_names))

pathway_names <- gsub("HALLMARK_", "", pathway_names)

# Select node names
node_names <- c("participant_id", "visit_month", "study", "diagnosis_at_baseline", "age_at_baseline", "sex", "race", "case_control_other_at_baseline", "global_famhistory", "GBA_mut", "LRRK2_mut"
                ,"SNCA_mut", "APOE_E4_mut", "PD_Mut", "levodopa", "dopamine", "other_med", "UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4", "MOCA", "Mobility39", "ADL39","Emotional39", "Stigma39", "Social39", "Cognition39","Communication39","Discomfort39", "Schwad_ADL", "REM", "ESS", "UPSIT"
                ,pathway_names)
colnames(df_no_controls) <- node_names

# Induce community structure, and create dataset which contains this information
community_structure <- c(rep("Demographics", 9),rep("Mutations", 5),rep("Medication", 3), rep("General PD Severity", 4), rep("Cognitive", 1) , rep("Disability", 9), rep("Sleep", 2), rep("Smell",1), rep("Gene expression",length(pathway_names)))
structure <- data.frame(node_names, community_structure)


###############################################################################
####################### Create stratified groups ##############################

dim(df_no_controls)
df_with_mutations <- df_no_controls %>% filter(PD_Mut == "Yes")
df_without_mutations <- df_no_controls %>% filter(PD_Mut != "Yes")
dim(df_with_mutations)
dim(df_without_mutations)

df_GBA <- df_no_controls %>% filter(GBA_mut == "Yes")
df_LRRK2 <- df_no_controls %>% filter(LRRK2_mut == "Yes")
df_SNCA <- df_no_controls %>% filter(SNCA_mut == "Yes")
df_APOE <- df_no_controls %>% filter(APOE_E4_mut == "Yes")

#2. Scale selected variables and order by ID

# Ordered by participant id
long_visit <- df_no_controls %>%
  arrange(participant_id) %>%
  group_by(participant_id)

# Scale to z-score continuous variables
# List of variables to scale
variables_to_scale <- c("UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4", "MOCA", "Mobility39", 
                        "ADL39", "Emotional39", "Stigma39", "Social39", "Cognition39", 
                        "Communication39", "Discomfort39", "Schwad_ADL", "ESS", "UPSIT",
                        pathway_names)

structure <- subset(structure, structure$node_names %in% variables_to_scale)

merged_df <- structure %>%
  left_join(NES, by = c("node_names" = "pathway"))  # Ajusta "pathway" al nombre correcto
merged_df$regulation[c(1:16)] <- merged_df$community_structure[c(1:16)]

# Loop through each variable and scale it
# for (variable in variables_to_scale) {
#   long_visit[[variable]] <- scale(long_visit[[variable]])
# }

#Check, patients are not filtered by idiopathic PD, all of them are included
table(long_visit$diagnosis_at_baseline) 


#3. Build multiple Linear Mixed Model 

#Select variables to include on network
dependent_variables <- variables_to_scale

filtered_variables <- structure$node_names[structure$community_structure != "Gene expression"]

variables_string <- paste(variables_to_scale, collapse = " + ")
variables_string_filtered <-  paste(filtered_variables, collapse = " + ")


library(lme4)
models <- list()
# Loop through each dependent variable
for (dependent_variable in dependent_variables) {
  # Create formula for the model
  #if (dependent_variable)
  
  if (structure$community_structure[structure$node_names==dependent_variable] == "Gene expression"){
    formula <- as.formula(paste(dependent_variable,  "~ ", variables_string_filtered, "+ (1|sex) + (1|participant_id) + (1|age_at_baseline)"))
  }
  else{
    formula <- as.formula(paste(dependent_variable, "~ ", variables_string, "+ (1|sex) + (1|participant_id) + (1|age_at_baseline) "))
  }
  # Create and fit the model
  model <- glmer(formula, data = long_visit)
  
  # Print the model summary
  #print(summary(model))
  
  # Optionally, you can store the model in a list or perform further analysis
  models[[dependent_variable]] <- summary(model)$coefficients[,3]
}

# Encontrar el número máximo de coeficientes entre todos los modelos
max_length <- max(sapply(models, length))

# Añadir ceros a los modelos que tienen menos coeficientes
for (depentent_variable in dependent_variables) {
  if (length(models[[depentent_variable]]) < max_length) {
    # Calcular cuántos ceros agregar
    diff <- max_length - length(models[[depentent_variable]])
    # Agregar los ceros
    models[[depentent_variable]] <- c(models[[depentent_variable]], rep(0, diff))
  }
}

# Convert list of lists into matrix
score<- do.call(rbind, models)
colnames(score)[1] <- variables_to_scale[1]

source("shifted_matrix.R")
# Change format of score matrix using external functions 
score_matrix <- score_matrix(score, shift_matrix(score))

# Plot adjusted node size
#qgraph(score_matrix, groups = structure$community_structure,  layout="spring", color=c("lightgreen", "lightblue","#00BFFF","orange", "pink","grey") ) 

# Addition ( M + t(M) )
symm_score <- (score_matrix + t(score_matrix) )/2

# Plot new network (symmetrical)
#qgraph(symm_score, groups = structure$community_structure,  layout="spring", color=c("lightgreen", "lightblue", "#00BFFF","orange", "pink","grey") ) 

spring_layout <- qgraph(symm_score, layout = "spring", DoNotPlot = TRUE)$layout

##Transform interaction between clinical variables to 0

is_gene <- structure$community_structure == "Gene expression"

symm_score_muted <- symm_score
# Bucle para convertir interacciones entre nodos que no son genes en 0 (eliminamos interacciones clinica-clinica)
for (i in 1:nrow(symm_score_muted)) {
  for (j in 1:ncol(symm_score_muted)) {
    # Si ambos nodos no son genes, establecer la interacción en 0
    if (!is_gene[i] & !is_gene[j]) {
      symm_score_muted[i, j] <- 0
    }
  }
}
# Plot new network (without gene-gene interaction)

qgraph(symm_score, groups = merged_df$regulation,  layout=spring_layout, color=c("lightgreen", "lightblue", "orange","coral","#00BFFF", "pink","grey"),
       labels = colnames(symm_score)) 

#png("~/Documents/Parkinson/GSEA/pathways_clinical.png", width = 2400, height = 1800, res = 300)  # 8x6 pulgadas a 300 dpi

qgraph(symm_score_muted, groups = merged_df$regulation,  layout=spring_layout, color=c("lightgreen", "lightblue", "orange","coral","#00BFFF", "pink","grey"),
       labels = colnames(symm_score_muted), label.cex = 1.3) 

dev.off()


setwd("~/Documents/Parkinson/transcriptomic//")

# Write output
#score_df <- as.data.frame(symm_score_muted)
#score_df<- cbind(row_names = rownames(score_df), score_df)
#write_csv(score_df, "transcriptomic_clinical_score_matrix.csv")

# Paso 1: Calcular la suma de interacciones por nodo (suma por fila o columna, es igual porque es simétrico)
node_strength <- rowSums(abs(symm_score_muted[pathway_names,]))

# Paso 2: Convertirlo a dataframe
strength_df <- data.frame(
  node = names(node_strength),
  strength = node_strength
)

# Paso 4: Filtrar solo "protein expression" y seleccionar los 10 más conectados
top10_genes <- rownames(strength_df[order(-strength_df$strength), ][1:10, ])


# Resultado: nodos de proteína más conectados
print(top10_genes)

selected_vars <- c("UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4", "MOCA", "Mobility39", 
                   "ADL39", "Emotional39", "Stigma39", "Social39", "Cognition39",
                   "Communication39", "Discomfort39", "Schwad_ADL", "ESS", "UPSIT",
                   top10_genes)

submatrix <- symm_score_muted[selected_vars, selected_vars]
structure <- subset(structure, structure$node_names %in% selected_vars)

merged_df <- structure %>%
  left_join(NES, by = c("node_names" = "pathway"))  # Ajusta "pathway" al nombre correcto

merged_df$regulation[c(1:16)] <- merged_df$community_structure[c(1:16)]

#To save image!!
setwd("~/Documents/Parkinson/GSEA")

# Open a PNG graphics device
png("pathways_clinical_top10.png", width = 1800, height = 1000, res = 500)

# Generate the qgraph plot
qgraph(submatrix,  vsize = 4,  groups = merged_df$regulation,  layout="spring", color=c("lightgreen", "lightblue", "orange","coral","#00BFFF", "pink","grey"),
       legend.cex = 0.2, labels = colnames(submatrix), label.cex = 1.5) 


# Close the graphics device
dev.off()

setwd("~/Documents/Parkinson/transcriptomic/")

score_output <- as.data.frame(submatrix)
score_output<- cbind(row_names = rownames(score_output), score_output)
write_csv(score_output, "transcriptomics_clinical_top10.csv")








################################################################################

#############################Mutations analysis#################################
#symm_scoreNOmutations <- symm_score_muted
symm_scoreNOmutations

#Mutations_score <- symm_score_muted
Mutations_score

library(MariNET)
difference <- differentiation(Mutations_score, symm_scoreNOmutations)
qgraph(symm_score_mutations, groups = structure$community_structure, layout=spring_layout, color=c("lightgreen", "lightblue", "#00BFFF","orange", "pink","grey"),
       labels = custom_labels) 


#saveRDS(symm_score_muted, "interactions_gene_clinical_muted.rds")

# Guardar la figura como PNG
png("Omics_integration/images_HD/network_comparison.png", width = 2000, height = 2000, res = 300)

layout_matrix <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE)
layout(layout_matrix)

# Adjust graphical parameters for reduced margins and larger titles
par(mar = c(0, 0, 0, 0))  # Set margins to zero
par(cex.main = 1)         # Increase title font size

# ------------------------- Plot 1: MariNET Network ------------------------------
qgraph(Mutations_score,
       title = "(A) Patients with any mutations",
       layout = spring_layout,
       groups = structure$community_structure,
       color = c("lightgreen", "lightblue", "#00BFFF","orange", "pink","grey"),
       labels = custom_labels,
       legend = FALSE,
       vsize = 7,
       label.font = 1.3)

# ----------------------------- Plot 2: Legend -----------------------------------
# Create an empty plot area for the legend
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
# Add a legend to the empty plot using unique community labels and colors
legend("center",
       legend = unique(structure$community_structure),
       fill =  c("orange","lightgreen", "lightblue", "pink","grey","#00BFFF"),
       title = "",
       cex = 1.1,    # Increase legend text size
       bty = "n")    # Remove the legend border

# ----------------------- Plot 3: EBICGlasso Network -------------------------------
qgraph(symm_scoreNOmutations,
       title = "(B) Patients with no mutations",
       layout = spring_layout,
       groups = structure$community_structure,
       color = c("lightgreen", "lightblue", "#00BFFF","orange", "pink","grey"),
       labels = custom_labels,
       legend = FALSE,
       vsize = 7,
       label.font = 1.3)

# ------------- Plot 4: Difference Between Methods Network -----------------------
qgraph(difference,
       title = "(C) Difference between both",
       layout = spring_layout,
       groups = structure$community_structure,
       color =  c("lightgreen", "lightblue", "#00BFFF","orange", "pink","grey"),
       labels = custom_labels,
       legend = FALSE,
       vsize = 7,
       label.font = 1.3)

# Reset the layout to the default single plot
layout(1)

# Finalizar guardado
dev.off()
