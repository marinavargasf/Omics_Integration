##################################################################################
#                                                                                #
#                        Parkinson's Disease Clinical Study                      #  
#                                                                                #
##################################################################################

# Set the working directory to the location of your project files
#setwd("~/Documents/Parkinson")

# --------------------- Run Required Previous Scripts --------------------------
# Source previous scripts that load data and perform preliminary processing.
# These scripts load clinical metadata, create a shifted matrix, and process data.
source("Data_processing/DataReading.R")
source("Data_processing/DataPreprocessing.R") 
# ------------------------- 1. Collect and Filter Data ---------------------------
# Merge clinical_data and metadata by participant_id using a left join. 
# This keeps only the clinical information for patients who have corresponding records.
master_matrix <- merge(clinical_data, metadata, by = "participant_id", all.x = TRUE)
dim(master_matrix)  # Check dimensions to see the number of merged records

# Reorder and select all variables used in the analysis including metadata information.
master_matrix <- master_matrix %>% select(
  participant_id, visit_month, study, diagnosis_at_baseline, age_at_baseline, sex, race, 
  case_control_other_at_baseline, global_famhistory,
  has_known_GBA_mutation_in_WGS, has_known_LRRK2_mutation_in_WGS,has_known_SNCA_mutation_in_WGS,
  has_known_APOE_E4_mutation_in_WGS, has_known_PD_Mutation_in_WGS,
  # test_abeta, test_ptau, test_tau, guid, study_participant_id  # Uncomment if needed
  mds_updrs_part_i_summary_score, mds_updrs_part_ii_summary_score, 
  mds_updrs_part_iii_summary_score, mds_updrs_part_iv_summary_score, 
  code_upd2hy_hoehn_and_yahr_stage, moca_total_score,
  pdq39_mobility_score, pdq39_adl_score, pdq39_emotional_score, pdq39_stigma_score, 
  pdq39_social_score, pdq39_cognition_score, pdq39_communication_score,
  pdq39_discomfort_score, Schwad_ADL_score, REM_score, ess_sleepiness_score, 
  upsit_total_score
)

#Selecting only patients!! with or without mutations
master_matrix <- master_matrix[master_matrix$diagnosis_at_baseline!="No PD Nor Other Neurological Disorder",]

# Filter by mutation
mut_cols <- grep("_mut$", colnames(master_matrix), value = TRUE)

# Replace NA values with "No" in the mutation columns
master_matrix[mut_cols] <- master_matrix[mut_cols] %>%
  lapply(function(x) replace(x, is.na(x), "No"))

# Filter rows where any "mut" column contains "Yes"
mutated_patients <- master_matrix %>%
  filter(rowSums(across(all_of(mut_cols), ~ . == "Yes")) > 0)

# Filter rows where "mut" column dont contain "Yes"
non_mutated_patients <- master_matrix %>%
  filter(rowSums(across(all_of(mut_cols), ~ . == "Yes")) == 0)

# ------------------------- Change Node Names ----------------------------------
# Define a vector of new variable names for clarity
node_names <- c("participant_id", "visit_month", "study", "diagnosis_at_baseline", 
                "age_at_baseline", "sex", "race", "case_control_other_at_baseline", 
                "global_famhistory","GBA_mutation", "LRRK2_mutation","SNCA_mutation",
                "APOE_E4_mutation", "has_known_PD_Mutation_in_WGS",
                "UPDRS1", "UPDRS2","UPDRS3", "UPDRS4", "HOEHN_Stage", "MOCA", "Mobility39", "ADL39",
                "Emotional39", "Stigma39", "Social39", "Cognition39", "Communication39",
                "Discomfort39", "Schwad_ADL", "REM", "ESS", "UPSIT")

# Rename the columns in master_matrix using the new node_names vector
colnames(master_matrix) <- node_names

# ----------------- Induce Community Structure for the Network ------------------
community_structure <- c(
  rep("Demographics", 9),
  rep("Genomics", 5),
  rep("General PD Severity", 5),
  rep("Cognitive", 1),
  rep("Disability", 9),
  rep("Sleep", 2),
  rep("Smell", 1)
)
structure <- data.frame(node_names, community_structure)

# --------------------- Filter by NA Values ------------------------------------
# Convert master_matrix to a data frame and filter out rows with more than 7 missing values
multiple_visit <- as.data.frame(master_matrix)
multiple_visit <- multiple_visit %>% filter(rowSums(is.na(.)) <= 7)

dim(multiple_visit)  # Check the dimensions after filtering (e.g., 4096 patients)

# Calculate the percentage of missing values for each variable
porcentaje_nas_variables <- as.data.frame(colMeans(is.na(multiple_visit)))

# Define a threshold (umbral) for allowed missing data (70% in this example)
umbral <- 0.7
# Select columns with a percentage of missing values below the threshold
columnas_filtradas <- colnames(multiple_visit)[porcentaje_nas_variables < umbral]
# Print the selected variable names
print(columnas_filtradas)

# Retain only the filtered variables in the final dataset
multiple_visit <- multiple_visit[, columnas_filtradas]

# Update the structure dataframe to keep only the filtered node names
structure <- subset(structure, structure$node_names %in% columnas_filtradas)

# ------------------- List Variables to Scale ----------------------------------
# Define which variables will be scaled for the subsequent network analyses
variables_to_scale <- c("UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4", "MOCA", "Mobility39",
                        "ADL39", "Emotional39", "Stigma39", "Social39", "Cognition39",
                        "Communication39", "Discomfort39", "Schwad_ADL", "ESS", "UPSIT")

# ----------------- Network Analysis with MariNET and EBICglasso ----------------
# Load MariNET package for network analysis based on linear mixed models
library("MariNET")

# Run a linear mixed model analysis (MariNET) on the filtered data
model <- lmm_analysis(multiple_visit, variables_to_scale)

# Estimate a network using the EBICglasso method with Spearman correlation
network <- estimateNetwork(multiple_visit[variables_to_scale],
                           default = "EBICglasso",
                           tuning = 0.25,
                           threshold = FALSE,
                           corMethod = "spearman")

# Compare the network derived from MariNET with the one from EBICglasso
difference_models <- differentiation(model, network$graph)

# ----------------- Extract Community Structure for Plotting ---------------------
# Obtain the community structure for the variables to scale from the structure dataset
structure_plot <- structure %>%
  filter(node_names %in% variables_to_scale) %>%
  pull(community_structure)

# -------------------- Prepare Custom Labels and Layout --------------------------
# Define custom labels for the network graph
custom_labels <- c("UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4", "MOCA", "Mob39",
                   "ADL39", "Emot39", "Stig39", "Social39", "Cogn39",
                   "Commu39", "Discomf39", "Schw_ADL", "ESS", "UPSIT")

# Load qgraph package (if not already loaded)
library(qgraph)

# Define a custom layout matrix for plotting multiple graphs:
# 1: Top-left (Graph 1)
# 2: Top-right (Graph 2)
# 3: Bottom-left (Graph 3)
# 4: Bottom-right (Legend)
# (Optional) To save high resolution png

png("genomics_integration.png", 
    width = 2400,       
    height = 2000,       
    res = 300,
    type = "cairo",
    bg = "transparent")  

layout_matrix <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE)
layout(layout_matrix)

# Adjust graphical parameters for reduced margins and larger titles
par(mar = c(0, 0, 0, 0))  # Set margins to zero
par(cex.main = 1)         # Increase title font size

# ------------------------- Plot 1: MariNET Network ------------------------------
edge_labels <- ifelse(abs(model) > 11, round(model, 2), NA)
edge_labels[5,14] = edge_labels[14,5] <- round(model[5,14],3) # to print specific interaction
edge_label_colors <- ifelse(edge_labels > 0, "#33CC33",
                            ifelse(edge_labels < 0, "red", "black"))

qgraph(model,
       title = "(A) MariNET",
       layout = "spring",
       groups = structure_plot,
       color = c("lightgreen", "lightblue", "orange", "pink", "grey"),
       edge.labels = edge_labels,     
       edge.label.cex = 1,   
       edge.label.font = 2,         # negrita
       edge.label.color = edge_label_colors,
       legend = FALSE, 
       nodeNames = names,
       vsize = 7, 
       label.cex = 1.2)      

# ----------------------------- Plot 2: Legend -----------------------------------
# Create an empty plot area for the legend
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
# Add a legend to the empty plot using unique community labels and colors
legend("center",
       legend = unique(structure_plot),
       fill = c("lightgreen", "lightblue", "orange", "pink", "grey"),
       title = "",
       cex = 1.1,    # Increase legend text size
       bty = "n")    # Remove the legend border

# ----------------------- Plot 3: EBICGlasso Network -------------------------------
edge_labels <- ifelse(abs(network$graph) > 0.3, round(network$graph, 2), NA)
edge_label_colors <- ifelse(edge_labels > 0, "#33CC33",
                            ifelse(edge_labels < 0, "red", "black"))

qgraph(network$graph,
       title = "(B) EBICGlasso",
       layout = "spring",
       groups = structure_plot,
       color = c("lightgreen", "lightblue", "orange", "pink", "grey"),
       edge.labels = edge_labels,     
       edge.label.cex = 1,   
       edge.label.font = 2,         # negrita
       edge.label.color = edge_label_colors,
       legend = FALSE, 
       nodeNames = names,
       vsize = 7, 
       label.cex = 1.2)      

# ------------- Plot 4: Difference Between Methods Network -----------------------
edge_labels <- ifelse(abs(difference_models) > 0.1, round(difference_models, 2), NA)
edge_label_colors <- ifelse(edge_labels > 0, "#33CC33",
                            ifelse(edge_labels < 0, "red", "black"))


qgraph(difference_models,
       title = "(C) Difference between methods",
       layout = "spring",
       groups = structure_plot,
       color = c("lightgreen", "lightblue", "orange", "pink", "grey"),
       edge.labels = edge_labels,     
       edge.label.cex = 1,   
       edge.label.font = 2,         # negrita
       edge.label.color = edge_label_colors,
       legend = FALSE, 
       nodeNames = names,
       vsize = 7, 
       label.cex = 1.2)      

# Reset the layout to the default single plot
layout(1)

dev.off()

# Comparison metrics 

W_mari   <- normalization(model)
W_glasso  <- network$graph
W_dif <- difference_models


metrics <- compare_three_networks(W_mari, W_glasso, W_dif)

print(metrics$PairwiseMetrics)

print(metrics$TopologicalMetrics)

