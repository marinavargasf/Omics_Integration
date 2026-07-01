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
data_normalized <- as.data.frame(readRDS("rnaseq/Results_mayo/Normalization/normalized_expression.rds")) #IMPORTANTE, la que esta guardada en Results_feb25 tiene visit 0.5
top_10_genes <- as.data.frame(readRDS("rnaseq/Results_mayo/Results_mayo/top10_genes.rds"))

rownames(data_normalized) <- sub("\\..*$", "", rownames(data_normalized))

# Load necessary library
library(biomaRt)

# Connect to Ensembl using biomaRt
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")  # Change dataset for other species

# Get Ensembl IDs from row names
ensembl_ids <- rownames(data_normalized)

# Map Ensembl IDs to Gene Names (HGNC symbols)
gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = mart
)

# Create a dataframe with Ensembl IDs and corresponding gene names
annotations <- data.frame(ensembl_gene_id = ensembl_ids, gene_name = NA)

# For each Ensembl ID in result_deseq, assign the corresponding gene name from gene_mapping
for (i in 1:nrow(annotations)) {
  gene_name <- gene_mapping$hgnc_symbol[gene_mapping$ensembl_gene_id == annotations$ensembl_gene_id[i]]
  
  # If a gene name is found, assign it; otherwise, leave as NA (Ensembl ID remains)
  if (length(gene_name) > 0) {
    annotations$gene_name[i] <- gene_name
  }
}

# Ensure that annotations has the same number of rows as result_deseq
if (nrow(annotations) != nrow(data_normalized)) {
  stop("Annotations and result_deseq have different numbers of rows")
}

# Now map the gene names to result_deseq (gene names from annotations)
data_normalized$gene_name <- annotations$gene_name

# Step 2: Handle cases where Gene Name is NA (replace with Ensembl ID)
data_normalized$gene_name[is.na(data_normalized$gene_name) | data_normalized$gene_name == ""] <- rownames(data_normalized)[is.na(data_normalized$gene_name) | data_normalized$gene_name == ""]

rownames(data_normalized) <- data_normalized$gene_name
data_normalized <- data_normalized[, !colnames(data_normalized) %in% c("gene_name")]

selected_genes <- rownames(top_10_genes)

df_filtered <- data_normalized[rownames(data_normalized) %in% selected_genes, ]

#Creation of gene_Expr matrix with sample_id, visit_month
gene_expr <- as.data.frame(t(df_filtered))

gene_expr$participant_id <- sub("^((.*?)-(.*?))-.*", "\\1", rownames(gene_expr))

gene_expr$time <- sub("^.*?-.*?-", "", rownames(gene_expr))

gene_expr$visit_month <- recode(gene_expr$time,
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



#colnames(data_normalized)[colnames(data_normalized) == "participant"] <- "participant_id"

# View(metadata)

clinical_data <- merge(clinical_data, gene_expr, by.x = c("participant_id", "visit_month"), by.y =c("participant_id", "visit_month"), all=TRUE)

master_matrix <- merge(clinical_data, metadata,by = "participant_id", all.x = T) 

# left join, only need clinical information on patients who have registers 
dim(master_matrix)

#Filter only cases to check influence on gene expression

df_no_controls <- master_matrix %>%
  filter(!tolower(trimws(case_control_other_at_baseline)) %in% c("control", "other")) #tolower pasa a minusculas

table(df_no_controls$case_control_other_at_baseline)

library(dplyr)

df_no_controls <- df_no_controls %>% select( participant_id, visit_month, study, diagnosis_at_baseline, age_at_baseline, sex, race, case_control_other_at_baseline, global_famhistory, has_known_GBA_mutation_in_WGS, has_known_LRRK2_mutation_in_WGS,
                                            has_known_SNCA_mutation_in_WGS, has_known_APOE_E4_mutation_in_WGS, has_known_PD_Mutation_in_WGS, on_levodopa, on_dopamine_agonist, on_other_pd_medications, mds_updrs_part_i_summary_score, mds_updrs_part_ii_summary_score,
                                            mds_updrs_part_iii_summary_score, mds_updrs_part_iv_summary_score, moca_total_score, pdq39_mobility_score, pdq39_adl_score, pdq39_emotional_score, pdq39_stigma_score, pdq39_social_score, pdq39_cognition_score, pdq39_communication_score, pdq39_discomfort_score, 
                                            Schwad_ADL_score, REM_score, ess_sleepiness_score, upsit_total_score, ANXA3, GPR27, FBXL13, LINC02656, HECW2, ACSL1, LINC01127, FOLR3, PLB1, MMP9)

# Select node names
node_names <- c("participant_id", "visit_month", "study", "diagnosis_at_baseline", "age_at_baseline", "sex", "race", "case_control_other_at_baseline", "global_famhistory", "GBA_mut", "LRRK2_mut"
                ,"SNCA_mut", "APOE_E4_mut", "PD_Mut", "levodopa", "dopamine", "other_med", "UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4", "MOCA", "Mobility39", "ADL39","Emotional39", "Stigma39", "Social39", "Cognition39","Communication39","Discomfort39", "Schwad_ADL", "REM", "ESS", "UPSIT"
                ,"ANXA3", "GPR27", "FBXL13", "LINC02656", "HECW2", "ACSL1", "LINC01127", "FOLR3", "PLB1", "MMP9")
colnames(df_no_controls) <- node_names

# Induce community structure, and create dataset which contains this information
community_structure <- c(rep("Demographics", 9),rep("Mutations", 5),rep("Medication", 3), rep("General PD Severity", 4), rep("Cognitive", 1) , rep("Disability", 9), rep("Sleep", 2), rep("Smell",1), rep("Gene expression",10))
structure <- data.frame(node_names, community_structure)

###############################################################################
####################### Calculate NA presence in data #########################

#Check mutations NA 
sum(is.na(mutations$has_known_GBA_mutation_in_WGS))
sum(is.na(mutations$has_known_PD_Mutation_in_WGS))



missing_summary <- df_no_controls %>%
  summarise(across(everything(), ~ mean(!is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "completeness")

library(ggplot2)

ggplot(missing_summary, aes(x = reorder(variable, completeness), y = completeness)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Porcentaje de presencia por variable",
       x = "Variable",
       y = "% de valores presentes (no NA)") +
  theme_minimal()


n_full_cases <- sum(complete.cases(df_no_controls))
cat("Casos completos actuales:", n_full_cases, "\n")

impact_table <- sapply(names(df_no_controls), function(var) {
  data_subset <- df_no_controls[, setdiff(names(df_no_controls), var)]
  sum(complete.cases(data_subset))
})

impact_df <- data.frame(
  variable = names(impact_table),
  full_cases_after_removal = as.integer(impact_table)
)

# Ordenar de mayor a menor impacto positivo
impact_df <- impact_df %>%
  arrange(desc(full_cases_after_removal))

print(impact_df)


###############################################################################
####################### Create stratified groups ##############################


dim(df_no_controls)
df_with_mutations <- df_no_controls %>% filter(PD_Mut == "Yes")
df_without_mutations <- df_no_controls %>% filter(PD_Mut == "No")
dim(df_with_mutations)
dim(df_without_mutations)

df_GBA <- df_no_controls %>% filter(GBA_mut == "Yes")
df_NO_GBA <- df_no_controls %>% filter(GBA_mut == "No")
#df_GBA_filtered <- na.omit((df_GBA))

df_LRRK2 <- df_no_controls %>% filter(LRRK2_mut == "Yes")
df_NO_LRRK2 <- df_no_controls %>% filter(LRRK2_mut == "No")
#df_LRRK2_filtered <- na.omit((df_LRRK2))

df_SNCA <- df_no_controls %>% filter(SNCA_mut == "Yes")
df_NO_SNCA <- df_no_controls %>% filter(SNCA_mut == "No")
#df_SNCA_filtered <- na.omit((df_SNCA))

df_APOE <- df_no_controls %>% filter(APOE_E4_mut == "Yes")
df_NO_APOE <- df_no_controls %>% filter(APOE_E4_mut == "No")
#df_APOE_filtered <- na.omit((df_APOE))


#2. Scale selected variables and order by ID

# Ordered by participant id
long_visit <- df_NO_GBA %>%
  arrange(participant_id) %>%
  group_by(participant_id)

# Scale to z-score continuous variables
# List of variables to scale
variables_to_scale <- c("UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4", "MOCA", "Mobility39", 
                        "ADL39", "Emotional39", "Stigma39", "Social39", "Cognition39", 
                        "Communication39", "Discomfort39", "Schwad_ADL", "ESS", "UPSIT",
                        "ANXA3", "GPR27", "FBXL13", "LINC02656", "HECW2", "ACSL1", "LINC01127", "FOLR3", "PLB1", "MMP9")

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

filtered_variables <- structure$node_names[structure$community_structure != "Gene expression"]

variables_string <- paste(variables_to_scale, collapse = " + ")
variables_string_filtered <-  paste(filtered_variables, collapse = " + ")


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
custom_labels <-  c("UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4", "MOCA", "Mob39", 
                    "ADL39", "Emot39", "Stig39", "Social39", "Cogn39", 
                    "Commu39", "Discomf39", "Schw_ADL", "ESS", "UPSIT",
                    "ANXA3", "GPR27", "FBXL13", "LINC02656", "HECW2", "ACSL1", "LINC01127", "FOLR3", "PLB1", "MMP9")

qgraph(symm_score, groups = structure$community_structure,  layout="spring", color=c("lightgreen", "lightblue", "#00BFFF","orange", "pink","grey"),
       labels = custom_labels) 

#png("Omics_integration/images_HD/genes_clinical_NOGBAmutations.png", width = 2400, height = 1800, res = 300)  # 8x6 pulgadas a 300 dpi

qgraph(symm_score_muted, groups = structure$community_structure, layout=spring_layout, color=c("lightgreen", "lightblue", "#00BFFF","orange", "pink","grey"),
       labels = custom_labels) 
#dev.off()

#GBA_mutations <- symm_score_muted
NO_GBA_mutations <- symm_score_muted
#APOE_mutations <- symm_score_muted


library(MariNET)
difference <- differentiation(GBA_mutations, NO_GBA_mutations)
qgraph(symm_score_mutations, groups = structure$community_structure, layout=spring_layout, color=c("lightgreen", "lightblue", "#00BFFF","orange", "pink","grey"),
       labels = custom_labels) 


#saveRDS(symm_score_muted, "interactions_gene_clinical_muted.rds")

# Guardar la figura como PNG
png("Omics_integration/images_HD/GBA_.png", width = 2000, height = 2000, res = 300)

layout_matrix <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE)
layout(layout_matrix)

# Adjust graphical parameters for reduced margins and larger titles
par(mar = c(0, 0, 0, 0))  # Set margins to zero
par(cex.main = 1)         # Increase title font size

# ------------------------- Plot 1: MariNET Network ------------------------------
qgraph(GBA_mutations,
       title = "(A) Patients with GBA mutation",
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
qgraph(NO_GBA_mutations,
       title = "(B) Patients with NO GBA mutation",
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
