################################################################################
# Automatized script for longitudinal analysis by linear mixed models
################################################################################


library(dplyr)      # filtering datasets
library(rstatix)    # summary statistics
library(ggpubr)     # convenient summary statistics and plots
library(lme4)       # linear mixed effects model
library(qgraph)     # network package
library(utils)

setwd(here::here("results/proteomics/"))

#Run required previous scripts
source(here::here("R/network_layout.R"))
source(here::here("clinical_processing/clinical_data.R"))


#1. Collect and filter data 
data_normalized <- as.data.frame(readRDS(here::here("results/proteomics/filtered_proteins_metadata_CSF_SYMBOL.rds"))) 
proteins <- colnames(data_normalized)[-c(1,2)]
proteins_1 <- proteins
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
                                          ,all_of(proteins))
# Select node names
node_names <- c("participant_id", "visit_month", "study", "diagnosis_at_baseline", "age_at_baseline",
                 "sex", "race", "case_control_other_at_baseline", "global_famhistory",
                 "levodopa", "dopamine", "other_med",
                 "GBA_mut", "LRRK2_mut", "SNCA_mut", "APOE_E4_mut", "PD_Mut",
                 "UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4", "MOCA", "Mobility39", "ADL39", "Emotional39",
                 "Stigma39", "Social39", "Cognition39", "Communication39", "Discomfort39", "Schwad_ADL",
                 "REM", "ESS", "UPSIT", proteins)
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
# for (variable in variables_to_scale) {
#   long_visit[[variable]] <- scale(long_visit[[variable]])
# }

#Check, patients are not filtered by idiopathic PD, all of them are included
table(long_visit$diagnosis_at_baseline) 


#3. Build multiple Linear Mixed Model (MariNET >= 2.0.0)

# Clinical variables are modeled on the rest of clinical variables and proteins,
# and proteins are modeled on clinical variables. Edge weight = mean of both t-values.
# Clinical-clinical and protein-protein interactions are set to 0 by the function.

clinical_variables <- setdiff(variables_to_scale, proteins)

symm_score <- MariNET::lmm_integration_analysis(
  as.data.frame(long_visit),               # Clinical data merged with CSF protein abundance
  variables_to_scale = clinical_variables, # Clinical variables (network nodes)
  molecular_features = proteins,              # Proteins (network nodes)
  random_effects = "(1|sex) + (1|participant_id) + (1|age_at_baseline)"
)

# Same node order as community structure (clinical variables, then proteins)
stopifnot(identical(colnames(symm_score), structure$node_names))

# Plot new network (without gene-gene interaction)
spring_layout <- network_layout(symm_score)

qgraph(symm_score,  groups = structure$community_structure,  layout=spring_layout, color=c("lightgreen", "lightblue","orange","#B9AEDC","pink","grey")  ) 

#To save image!!
setwd(here::here("results/proteomics/"))

# ------------------------------------------------------------------------
# CSF_proteins_clinical plot -- PNG, SVG, PDF
# ------------------------------------------------------------------------
png("CSF_proteins_clinical.png", width = 2400, height = 1800, res = 500)
qgraph(symm_score, vsize = 4,  groups = structure$community_structure,  layout=spring_layout, color=c("lightgreen", "lightblue","orange","#B9AEDC","pink","grey"), legend.cex = 0.2, labels = colnames(symm_score)) 
dev.off()

svg("CSF_proteins_clinical.svg", width = 2400, height = 1800)
qgraph(symm_score, vsize = 4,  groups = structure$community_structure,  layout=spring_layout, color=c("lightgreen", "lightblue","orange","#B9AEDC","pink","grey"), legend.cex = 0.2, labels = colnames(symm_score)) 
dev.off()

pdf("CSF_proteins_clinical.pdf", width = 2400, height = 1800)
qgraph(symm_score, vsize = 4,  groups = structure$community_structure,  layout=spring_layout, color=c("lightgreen", "lightblue","orange","#B9AEDC","pink","grey"), legend.cex = 0.2, labels = colnames(symm_score)) 
dev.off()

# Write output
score_df <- as.data.frame(symm_score)
score_df<- cbind(row_names = rownames(score_df), score_df)
write.csv(score_df, "CSF_proteins_clinical_score_matrix.csv")

###############################################################################
##############Para graficar solamente los top10 interactores####################
################################################################################


# Paso 1: Calcular la suma de interacciones por nodo (suma por fila o columna, es igual porque es simétrico)
node_strength <- rowSums(abs(symm_score[proteins,]))
node_max_strength <- apply(abs(symm_score[proteins, ]), 1, max)

# Paso 2: Convertirlo a dataframe
#strength_df <- data.frame(
#  node = names(node_strength),
#  strength = node_strength
#)
strength_df <- data.frame(
  node = names(node_max_strength),
  strength = node_max_strength
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


#To save image!!
setwd(here::here("results/proteomics/"))

# ------------------------------------------------------------------------
# CSF_proteins_clinical_top10 plot -- PNG, SVG, PDF (same physical size)
# ------------------------------------------------------------------------

el <- make_edge_labels(submatrix, threshold = 1.6)
ec_labels <- make_edge_colors(el)
ec_lines  <- make_edge_colors_gradient(submatrix, exponent = 1.6)

fig_w <- 10     # inches
fig_h <- 7.5    # inches

plot_net <- function(mat, layout, groups = structure$community_structure) {
  qgraph(mat,
         vsize = 4,
         groups = groups,
         layout = layout,
         color = c("lightgreen", "lightblue", "orange", "#B9AEDC", "pink", "grey"),
         legend.cex = 0.45,
         labels = colnames(mat),
         edge.labels = el,
         edge.label.cex = 0.7,
         edge.label.font = 2,
         edge.label.color = ec_labels,
         edge.color = ec_lines,
         fade = FALSE,
         label.cex = 1.2)
}

png("CSF_proteins_clinical_top10.png", width = fig_w, height = fig_h, units = "in", res = 300)
plot_net(submatrix, network_layout(submatrix)); dev.off()


svg("CSF_proteins_clinical_top10.svg", width = fig_w, height = fig_h)
plot_net(submatrix, network_layout(submatrix)); dev.off()


pdf("CSF_proteins_clinical_top10.pdf", width = fig_w, height = fig_h)
plot_net(submatrix, network_layout(submatrix)); dev.off()


setwd(here::here("results/proteomics/"))

score_output <- as.data.frame(submatrix)
score_output<- cbind(row_names = rownames(score_output), score_output)
write.csv(score_output, "CSF_proteins_clinical_top10.csv")

################################################################################
##############Plot only common proteins with CSF################################
################################################################################
library(readr)

setwd(here::here("results/proteomics/"))
# Step 1: Load CSF score matrix to extract CSF protein names
pla_normalized <- as.data.frame(readRDS(here::here("results/proteomics/filtered_proteins_metadata_PLA_SYMBOL.rds"))) 
#csf_proteins<- colnames(data_normalized)[-c(1,2)]

# top_prots_1<-readRDS("proteomic/differential_expression_proteins.rds")

# Plasma DAPs taken from the plasma MariNET input (output of 09), so this
# script does not depend on the plasma network script and runs first
pla_proteins <- colnames(pla_normalized)[-c(1, 2)]

# Step 2: Find ALL common proteins between plasma and CSF
common_proteins <- intersect(proteins, pla_proteins)
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

structure <- data.frame(node_names, community_structure)

structure <- subset(structure, structure$node_names %in% variables_to_scale)

# Loop through each variable and scale it
#for (variable in variables_to_scale) {
#  long_visit[[variable]] <- scale(long_visit[[variable]])
#}

#Check, patients are not filtered by idiopathic PD, all of them are included
table(long_visit$diagnosis_at_baseline) 


#3. Build multiple Linear Mixed Model (MariNET >= 2.0.0)

# Clinical variables are modeled on the rest of clinical variables and proteins,
# and proteins are modeled on clinical variables. Edge weight = mean of both t-values.
# Clinical-clinical and protein-protein interactions are set to 0 by the function.

clinical_variables <- setdiff(variables_to_scale, common_proteins)

symm_score <- MariNET::lmm_integration_analysis(
  as.data.frame(long_visit),               # Clinical data merged with CSF protein abundance
  variables_to_scale = clinical_variables, # Clinical variables (network nodes)
  molecular_features = common_proteins,       # Proteins (network nodes)
  random_effects = "(1|sex) + (1|participant_id) + (1|age_at_baseline)"
)

# Same node order as community structure (clinical variables, then proteins)
stopifnot(identical(colnames(symm_score), structure$node_names))

# Plot new network (without gene-gene interaction)
qgraph(symm_score,  groups = structure$community_structure,  layout="spring", color=c("lightgreen", "lightblue","orange","#B9AEDC","pink","grey")  ) 

# Step 7: Plot
setwd(here::here("results/proteomics/"))

el <- make_edge_labels(symm_score, threshold = 2)
ec_labels <- make_edge_colors(el)                     # for label text (unchanged)
ec_lines  <- make_edge_colors_gradient(symm_score, exponent = 1.6)      # for edge lines (independent

# ------------------------------------------------------------------------
# PLA_commonCSF_allproteins plot -- PNG, SVG, PDF
# ------------------------------------------------------------------------
lay <- network_layout(symm_score)

png("CSF_commonPLA_allproteins.png", width = fig_w, height = fig_h, units = "in", res = 300)
plot_net(symm_score, lay); dev.off()

svg("CSF_commonPLA_allproteins.svg", width = fig_w, height = fig_h)
plot_net(symm_score, lay); dev.off()

pdf("CSF_commonPLA_allproteins.pdf", width = fig_w, height = fig_h)
plot_net(symm_score, lay); dev.off()


# Step 8: Save matrix
score_output_common <- as.data.frame(symm_score)
score_output_common <- cbind(row_names = rownames(score_output_common), score_output_common)
write.csv(score_output_common, "CSF_commonPLA_allproteins_score_matrix.csv")