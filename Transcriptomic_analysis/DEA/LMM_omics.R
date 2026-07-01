################################################################################

#Automatized script for longitudinal analysis by linear mixed models

#Marina Vargas-Fernández


library(dplyr)      # filtering datasets
library(rstatix)    # summary statistics
library(ggpubr)     # convenient summary statistics and plots
library(tidyverse)  # not sure if usefull
library(lme4)       # linear mixed effects model
library(qgraph)     # network package

setwd("~/Documents/Parkinson")

#Run required previous scripts

source("clinical_metadata_v4_2023/DataReading_v4.R")
source("shifted_matrix.R")
source("NetworkComparison_3.R") #revisar este
source("Long_analysis_v4.R")

#1. Collect and filter data 

# View(metadata)

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
                                          ,pdq39_discomfort_score, Schwad_ADL_score, REM_score, ess_sleepiness_score, upsit_total_score)
# Select node names
node_names <- c("participant_id", "visit_month", "study", "diagnosis_at_baseline", "age_at_baseline", "sex", "race", "case_control_other_at_baseline", "global_famhistory", "GBA_mut", "LRRK2_mut"
                ,"SNCA_mut", "APOE_E4_mut", "PD_Mut", "levodopa", "dopamine", "other_med", "UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4", "MOCA", "Mobility39", "ADL39","Emotional39", "Stigma39", "Social39", "Cognition39","Communication39","Discomfort39", "Schwad_ADL", "REM", "ESS", "UPSIT")
colnames(master_matrix) <- node_names

# Induce community structure, and create dataset which contains this information
community_structure <- c(rep("Demographics", 9),rep("Mutations", 5),rep("Medication", 3), rep("General PD Severity", 4), rep("Cognitive", 1) , rep("Disability", 9), rep("Sleep", 2), rep("Smell",1) )
structure <- data.frame(node_names, community_structure)



#test for unfiltered: 

multiple_visit<-master_matrix
dim(multiple_visit) #27230 patients


#2. Scale selected variables and order by ID

# Ordered by participant id
long_visit <- multiple_visit %>%
  arrange(participant_id) %>%
  group_by(participant_id)

summary(long_visit)

#2.1 Modification of dataset to extreme values

long_visit <- long_visit %>%
  mutate(UPSIT = ifelse(sex == "Male" & UPDRS3>10, 40 , UPSIT))

#2.2 Visualization of modified variables 

  ggplot(long_visit, aes(x = sex, y = UPSIT, fill = sex)) +
    geom_boxplot() +
    labs(x = "Sex", y = "UPSIT") +
    ggtitle("Boxplot of UPSIT Score by Sex")
  
  ggplot(data = long_visit, aes(x = UPSIT, y = UPDRS3)) +
    geom_point()
  
long_visit_male <- filter(long_visit, sex == "Male")
  ggplot(data = long_visit_male, aes(x = UPSIT, y = UPDRS3)) +
    geom_point()
  

# Scale to z-score continuous variables
# List of variables to scale
variables_to_scale <- c("UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4", "MOCA", "Mobility39", 
                        "ADL39", "Emotional39", "Stigma39", "Social39", "Cognition39", 
                        "Communication39", "Discomfort39", "Schwad_ADL", "ESS", "UPSIT")

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
models <- list()
# Loop through each dependent variable
for (dependent_variable in dependent_variables) {
  # Create formula for the model
  formula <- as.formula(paste(dependent_variable, "~ UPDRS1 + UPDRS2 + UPDRS3 + UPDRS4 + MOCA + Mobility39 + ADL39 + Emotional39 + Stigma39 + Social39 + Cognition39 + Communication39 + Discomfort39 + Schwad_ADL + ESS + UPSIT + (1|participant_id) + (1|sex)"))
  
  # Create and fit the model
  model <- lmer(formula, data = long_visit)
  
  # Print the model summary
  #print(summary(model))
  
  # Optionally, you can store the model in a list or perform further analysis
  models[[dependent_variable]] <- summary(model)$coefficients[,3]
}

# Convert list of lists into matrix
score<- do.call(rbind, models)
colnames(score)[1] <- variables_to_scale[1]

# Change format of score matrix using external functions 
source("shifted_matrix.R")
score_matrix <- score_matrix(score, shift_matrix(score))

# Get node degree
grado <- rowSums(score_matrix) + colSums(score_matrix)

# Normalize degree
grado_normalizado <- grado / max(grado)

node_size <- 5 + 3.5 * grado_normalizado  # Adjust normalization to size

# Plot adjusted node size
#qgraph(score_matrix, vsize = node_size, groups = structure$community_structure,  layout="spring", color=c("lightgreen", "lightblue","orange", "pink","grey") ) 

# Addition ( M + t(M) )
symm_score <- score_matrix + t(score_matrix)

# Plot new network (symmetrical)
qgraph(symm_score, vsize = node_size, groups = structure$community_structure,  layout="spring", color=c("lightgreen", "lightblue","orange", "pink","grey") ) 


# Write output
#score_df <- as.data.frame(symm_score)
#write_csv(score_df, "score_matrix.csv")

#TEST TO SEE T-VALUE OUTPUT ON LME4 MODEL 
#formula <- "UPDRS1 ~ UPDRS2 + UPDRS3 + UPDRS4 + MOCA + Mobility39 + ADL39 + Emotional39 + Stigma39 + Social39 + Cognition39 + Communication39 + Discomfort39 + Schwad_ADL + ESS + UPSIT + (1|participant_id) + (1|sex)"
#summary(model)$coefficients


###############################################################################################
################Test to see if network estimation produces same result#########################
library(bootnet)

network <- estimateNetwork(long_visit[, c(18:31,33,34)], default = "EBICglasso", tuning = 0.25, threshold = FALSE, corMethod ="spearman")
qgraph(network$graph, groups = structure$community_structure,  layout="spring", color=c("lightgreen", "lightblue","orange", "pink","grey") ) 

################################################################################
##########################OMICS DATA############################################
################################################################################

omics <- readRDS("marina_trial.rds")
