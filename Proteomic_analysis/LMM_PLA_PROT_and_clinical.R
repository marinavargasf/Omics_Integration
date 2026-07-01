################################################################################
# Automatized script for longitudinal analysis by linear mixed models
################################################################################


library(dplyr)      # filtering datasets
library(rstatix)    # summary statistics
library(ggpubr)     # convenient summary statistics and plots
library(lme4)       # linear mixed effects model
library(qgraph)     # network package
library(utils)

setwd("~/Documents/Omics_Integration/Proteomic_analysis/")

#Run required previous scripts

source("~/Documents/Parkinson/clinical_metadata_v4_2023/DataReading_v4.R")
source("~/Documents/Parkinson/data_processing.R") 


#1. Collect and filter data 
data_normalized <- as.data.frame(readRDS("~/Documents/Omics_Integration/Proteomic_analysis/Results/filtered_proteins_metadata_PLA_SYMBOL.rds")) 
proteins <- colnames(data_normalized)[-c(1,2)]

# View(metadata)

clinical_data <- merge(clinical_data, data_normalized, by.x = c("participant_id", "visit_month"), by.y =c("participant_id", "visit_month"), all=TRUE)

master_matrix <- merge(clinical_data, metadata,by = "participant_id", all.x = T) 

# left join, only need clinical information on patients who have registers 
dim(master_matrix)
# All variables used in the analysis, including metadata information
master_matrix <- master_matrix %>% select(participant_id, visit_month, study, diagnosis_at_baseline, age_at_baseline, sex, race, case_control_other_at_baseline, global_famhistory
                                          ,on_levodopa, on_dopamine_agonist, on_other_pd_medications, has_known_GBA_mutation_in_WGS, has_known_LRRK2_mutation_in_WGS
                                          ,has_known_SNCA_mutation_in_WGS, has_known_APOE_E4_mutation_in_WGS, has_known_PD_Mutation_in_WGS
                                          #, test_abeta, test_ptau, test_tau, guid, study_participant_id
                                          ,mds_updrs_part_i_summary_score, mds_updrs_part_ii_summary_score, mds_updrs_part_iii_summary_score, mds_updrs_part_iv_summary_score, moca_total_score
                                          ,pdq39_mobility_score, pdq39_adl_score, pdq39_emotional_score, pdq39_stigma_score, pdq39_social_score, pdq39_cognition_score, pdq39_communication_score
                                          ,pdq39_discomfort_score, Schwad_ADL_score, REM_score, ess_sleepiness_score, upsit_total_score
                                          ,c(proteins))
# Select node names
node_names <- c("participant_id", "visit_month", "study", "diagnosis_at_baseline", "age_at_baseline", "sex", "race", "case_control_other_at_baseline", "global_famhistory", "GBA_mut", "LRRK2_mut"
                ,"SNCA_mut", "APOE_E4_mut", "PD_Mut", "levodopa", "dopamine", "other_med", "UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4", "MOCA", "Mobility39", "ADL39","Emotional39", "Stigma39", "Social39", "Cognition39","Communication39","Discomfort39", "Schwad_ADL", "REM", "ESS", "UPSIT"
                ,proteins)
colnames(master_matrix) <- node_names

# Induce community structure, and create dataset which contains this information
community_structure <- c(rep("Demographics", 9),rep("Mutations", 5),rep("Medication", 3), rep("General PD Severity", 4), rep("Cognitive", 1) , rep("Disability", 9), rep("Sleep", 2), rep("Smell",1), rep("Protein expression", length(proteins)))
structure <- data.frame(node_names, community_structure)



#test for unfiltered: 

multiple_visit<-master_matrix
dim(multiple_visit) #27397 patients

#1.2 Filter only diagnosis PD

table(multiple_visit$diagnosis_at_baseline)

multiple_visit <- subset(multiple_visit, diagnosis_at_baseline =="Idiopathic PD")

#2. Scale selected variables and order by ID

# Ordered by participant id
long_visit <- multiple_visit %>%
  arrange(participant_id) %>%
  group_by(participant_id)

# Scale to z-score continuous variables
# List of variables to scale
variables_to_scale <- c("UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4", "MOCA", "Mobility39", 
                        "ADL39", "Emotional39", "Stigma39", "Social39", "Cognition39", 
                        "Communication39", "Discomfort39", "Schwad_ADL", "ESS", "UPSIT",
                        proteins)

structure <- subset(structure, structure$node_names %in% variables_to_scale)


# Loop through each variable and scale it
for (variable in variables_to_scale) {
  long_visit[[variable]] <- scale(long_visit[[variable]])
}

#Check, patients are not filtered by idiopathic PD, all of them are included
table(long_visit$diagnosis_at_baseline) 


#3. Build multiple Linear Mixed Model 

#Select variables to include on network
dependent_variables <- variables_to_scale

filtered_variables <- structure$node_names[structure$community_structure != "Protein expression"]

variables_string <- paste(variables_to_scale, collapse = " + ")
variables_string_filtered <-  paste(filtered_variables, collapse = " + ")


models <- list()
# Loop through each dependent variable
for (dependent_variable in dependent_variables) {
  # Create formula for the model
  #if (dependent_variable)
  
  if (structure$community_structure[structure$node_names==dependent_variable] == "Protein expression"){
    formula <- as.formula(paste(dependent_variable,  "~ ", variables_string_filtered, "+ (1|participant_id)"))
  }
  else{
    formula <- as.formula(paste(dependent_variable, "~ ", variables_string, "+ (1|participant_id)"))
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

source("~/Documents/Parkinson/shifted_matrix.R")

# Change format of score matrix using external functions 
score_matrix <- score_matrix(score, shift_matrix(score))


# Plot adjusted node size
qgraph(score_matrix, groups = structure$community_structure,  layout="spring", color=c("lightgreen", "lightblue","orange","#B9AEDC", "pink","grey") ) 

# Addition ( M + t(M) )
symm_score <- (score_matrix + t(score_matrix) )/2

# Plot new network (symmetrical)
qgraph(symm_score, groups = structure$community_structure,  layout="spring", color=c("lightgreen", "lightblue","orange","#B9AEDC","pink","grey") ) 

##Transform interaction between clinical variables to 0

is_gene <- structure$community_structure == "Protein expression"

# Bucle para convertir interacciones entre nodos que no son genes en 0 (eliminamos interacciones clinica-clinica)
for (i in 1:nrow(symm_score)) {
  for (j in 1:ncol(symm_score)) {
    # Si ambos nodos no son genes, establecer la interacción en 0
    if (!is_gene[i] & !is_gene[j]) {
      symm_score[i, j] <- 0
    }
  }
}
# Plot new network (without gene-gene interaction)
qgraph(symm_score,  groups = structure$community_structure,  layout="spring", color=c("lightgreen", "lightblue","orange","#B9AEDC","pink","grey")  ) 


#To save image!!
setwd("~/Documents/Omics_Integration/Proteomic_analysis/Results/")

# Open a PNG graphics device
png("PLA_proteins_clinical.png", width = 1800, height = 1000, res = 500)

# Generate the qgraph plot
qgraph(symm_score, vsize = 3,  groups = structure$community_structure,  layout="spring", color=c("lightgreen", "lightblue","orange","#B9AEDC","pink","grey"), legend.cex = 0.2, labels = colnames(symm_score)) 

# Close the graphics device
dev.off()

# Write output
score_df <- as.data.frame(symm_score)
score_df<- cbind(row_names = rownames(score_df), score_df)
write.csv(score_df, "PLA_proteins_clinical_score_matrix.csv")


################################################################################
##############Para graficar solamente los top10 interactores####################
################################################################################


# Paso 1: Calcular la suma de interacciones por nodo (suma por fila o columna, es igual porque es simétrico)
node_strength <- rowSums(abs(symm_score[proteins,]))

# Paso 2: Convertirlo a dataframe
strength_df <- data.frame(
  node = names(node_strength),
  strength = node_strength
)

# Paso 4: Filtrar solo "protein expression" y seleccionar los 10 más conectados
top10_proteins <- rownames(strength_df[order(-strength_df$strength), ][1:10, ])


# Resultado: nodos de proteína más conectados
print(top10_proteins)

selected_vars <- c("UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4", "MOCA", "Mobility39", 
                   "ADL39", "Emotional39", "Stigma39", "Social39", "Cognition39",
                   "Communication39", "Discomfort39", "Schwad_ADL", "ESS", "UPSIT",
                   top10_proteins)

submatrix <- symm_score[selected_vars, selected_vars]
structure <- subset(structure, structure$node_names %in% selected_vars)

# Open a PNG graphics device
png("PLA_proteins_clinical_top10.png", width = 1800, height = 1000, res = 500)

# Generate the qgraph plot
qgraph(submatrix,  vsize = 4,  groups = structure$community_structure,  layout=layout_fixed, color=c("lightgreen", "lightblue","orange","#B9AEDC","pink","grey"), legend.cex = 0.2, labels = colnames(submatrix)) 


# Close the graphics device
dev.off()

score_output <- as.data.frame(submatrix)
score_output<- cbind(row_names = rownames(score_output), score_output)
write.csv(score_output, "PLA_proteins_clinical_top10.csv")

# score_output <- read_csv("PLA_proteins_clinical_top10.csv")

################################################################################
##############Plot only common proteins with CSF################################
################################################################################
library(readr)

setwd("~/Documents/Omics_Integration/Proteomic_analysis/")
# Step 1: Load CSF score matrix to extract CSF protein names
csf_normalized <- as.data.frame(readRDS("~/Documents/Omics_Integration/Proteomic_analysis/Results/filtered_proteins_metadata_CSF_SYMBOL.rds")) 
csf_proteins<- colnames(data_normalized)[-c(1,2)]

# top_prots_1<-readRDS("proteomic/differential_expression_proteins.rds")

csf_score <- read.csv("~/Documents/Omics_Integration/Proteomic_analysis/Results/CSF_proteins_clinical_score_matrix.csv", row.names=1)
csf_score <- csf_score[, -1] #Duplicated row 
csf_proteins <- colnames(csf_score)[-(1:16)]  # remove clinical vars

# Step 2: Find ALL common proteins between plasma and CSF
common_proteins <- intersect(proteins, csf_proteins)
print(paste("Common proteins found:", length(common_proteins)))
print(common_proteins)

# filtered_proteins_metadata_PLA <- readRDS("Results/filtered_proteins_metadata_PLA.rds")
# filtered_proteins_metadata_CSF <- readRDS("Results/filtered_proteins_metadata_CSF.rds")
# 
# common_proteins_HGNC <- intersect(colnames(filtered_proteins_metadata_PLA),colnames(filtered_proteins_metadata_CSF))
# 
# 
# 
# de_proteins_metadata_PLA <- readRDS("Results/differential_expression_proteins_PLA.rds")
# de_proteins_metadata_CSF <- readRDS("Results/differential_expression_proteins_CSF.rds")
# 
# common_proteins_DEG <- intersect(rownames(de_proteins_metadata_PLA),rownames(de_proteins_metadata_CSF))

# Step 3: Build selected variables with ALL common proteins

# Ordered by participant id
long_visit <- multiple_visit %>%
  arrange(participant_id) %>%
  group_by(participant_id)

# Scale to z-score continuous variables
# List of variables to scale
variables_to_scale <- c("UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4", "MOCA", "Mobility39", 
                        "ADL39", "Emotional39", "Stigma39", "Social39", "Cognition39", 
                        "Communication39", "Discomfort39", "Schwad_ADL", "ESS", "UPSIT",
                        common_proteins)

structure <- subset(structure, structure$node_names %in% variables_to_scale)

# Loop through each variable and scale it
for (variable in variables_to_scale) {
  long_visit[[variable]] <- scale(long_visit[[variable]])
}

#Check, patients are not filtered by idiopathic PD, all of them are included
table(long_visit$diagnosis_at_baseline) 


#3. Build multiple Linear Mixed Model 

#Select variables to include on network
dependent_variables <- variables_to_scale

filtered_variables <- structure$node_names[structure$community_structure != "Protein expression"]

variables_string <- paste(variables_to_scale, collapse = " + ")
variables_string_filtered <-  paste(filtered_variables, collapse = " + ")


models <- list()
# Loop through each dependent variable
for (dependent_variable in dependent_variables) {
  # Create formula for the model
  #if (dependent_variable)
  
  if (structure$community_structure[structure$node_names==dependent_variable] == "Protein expression"){
    formula <- as.formula(paste(dependent_variable,  "~ ", variables_string_filtered, "+ (1|participant_id)"))
  }
  else{
    formula <- as.formula(paste(dependent_variable, "~ ", variables_string, "+ (1|participant_id)"))
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

source("~/Documents/Parkinson/shifted_matrix.R")

# Change format of score matrix using external functions 
score_matrix <- score_matrix(score, shift_matrix(score))


# Plot adjusted node size
qgraph(score_matrix, groups = structure$community_structure,  layout="spring", color=c("lightgreen", "lightblue","orange","#B9AEDC", "pink","grey") ) 

# Addition ( M + t(M) )
symm_score <- (score_matrix + t(score_matrix) )/2

# Plot new network (symmetrical)
qgraph(symm_score, groups = structure$community_structure,  layout="spring", color=c("lightgreen", "lightblue","orange","#B9AEDC","pink","grey") ) 

##Transform interaction between clinical variables to 0

is_gene <- structure$community_structure == "Protein expression"

# Bucle para convertir interacciones entre nodos que no son genes en 0 (eliminamos interacciones clinica-clinica)
for (i in 1:nrow(symm_score)) {
  for (j in 1:ncol(symm_score)) {
    # Si ambos nodos no son genes, establecer la interacción en 0
    if (!is_gene[i] & !is_gene[j]) {
      symm_score[i, j] <- 0
    }
  }
}
# Plot new network (without gene-gene interaction)
qgraph(symm_score,  groups = structure$community_structure,  layout="spring", color=c("lightgreen", "lightblue","orange","#B9AEDC","pink","grey")  ) 

# Step 7: Plot
setwd("~/Documents/Omics_Integration/Proteomic_analysis/Results/")

png("PLA_commonCSF_allproteins.png", width = 1800, height = 1200, res = 500)
qgraph(symm_score,
       vsize = 4,
       groups = structure_common_full$community_structure,
       layout = layout_fixed_2,
       color = c("lightgreen", "lightblue", "orange", "#B9AEDC", "pink", "grey"),
       legend.cex = 0.2,
       labels = colnames(symm_score))
dev.off()

# Step 8: Save matrix
score_output_common <- as.data.frame(symm_score)
score_output_common <- cbind(row_names = rownames(score_output_common), score_output_common)
write.csv(score_output_common, "PLA_commonCSF_allproteins_score_matrix.csv")
