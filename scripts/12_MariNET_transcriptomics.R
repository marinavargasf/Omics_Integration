################################################################################

#Automatized script for longitudinal analysis by linear mixed models

#Marina Vargas-Fernández


library(dplyr)      # filtering datasets
library(rstatix)    # summary statistics
library(ggpubr)     # convenient summary statistics and plots
library(lme4)       # linear mixed effects model
library(qgraph)     # network package
library(MariNET)
stopifnot(packageVersion("MariNET") >= "2.0.0")

#Run required previous scripts

source(here::here("clinical_processing/clinical_data.R"))
source(here::here("R/network_layout.R"))

#1. Collect and filter data 
data_normalized <- as.data.frame(readRDS(here::here("results/transcriptomics/deseq2/normalized_expression_genesymbol.rds")))

dif_exp <- as.data.frame(readRDS(here::here("results/transcriptomics/deseq2/dif_exp_genesymbol_sig.rds")))

selected_genes <- rownames(dif_exp)

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

#saveRDS(gene_expr, "transcriptomic/gene_expr_longformat.rds")

colnames(gene_expr)[colnames(gene_expr) == "H1-0"] <- "H10"
selected_genes[selected_genes == "H1-0"] <- "H10"
#colnames(data_normalized)[colnames(data_normalized) == "participant"] <- "participant_id"

#To check how many samples we have
samples_full <- merge(gene_expr, metadata, by ="participant_id", all.x = T)
# To check how many cases and how many controls
samples_filtrado <- samples_full %>%
  group_by(participant_id) %>%
  slice(1) %>%
  ungroup()

table(samples_filtrado$diagnosis_at_baseline)
# View(metadata)

clinical_data_2 <- merge(clinical_data, gene_expr, by.x = c("participant_id", "visit_month"), by.y =c("participant_id", "visit_month"), all=TRUE)

master_matrix <- merge(clinical_data_2, metadata,by = "participant_id", all.x = T) 

# left join, only need clinical information on patients who have registers 
dim(master_matrix)

#Filter only cases to check influence on gene expression

df_no_controls <- subset(master_matrix, diagnosis_at_baseline =="Idiopathic PD")


table(df_no_controls$case_control_other_at_baseline)

library(dplyr)

df_no_controls <- df_no_controls %>% dplyr::select( participant_id, visit_month, study, diagnosis_at_baseline, age_at_baseline, sex, race, case_control_other_at_baseline, global_famhistory,
                                            on_levodopa, on_dopamine_agonist, on_other_pd_medications, mds_updrs_part_i_summary_score, mds_updrs_part_ii_summary_score,
                                            mds_updrs_part_iii_summary_score, mds_updrs_part_iv_summary_score, moca_total_score, pdq39_mobility_score, pdq39_adl_score, pdq39_emotional_score, pdq39_stigma_score, pdq39_social_score, pdq39_cognition_score, pdq39_communication_score, pdq39_discomfort_score, 
                                            Schwad_ADL_score, REM_score, ess_sleepiness_score, upsit_total_score,
                                            all_of(selected_genes))
# Select node names
node_names <- c("participant_id", "visit_month", "study", "diagnosis_at_baseline", "age_at_baseline", "sex", "race", "case_control_other_at_baseline", "global_famhistory"
                , "levodopa", "dopamine", "other_med", "UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4", "MOCA", "Mobility39", "ADL39","Emotional39", "Stigma39", "Social39", "Cognition39","Communication39","Discomfort39", "Schwad_ADL", "REM", "ESS", "UPSIT"
                ,selected_genes)
colnames(df_no_controls) <- node_names

# Induce community structure, and create dataset which contains this information
community_structure <- c(rep("Demographics", 9),rep("Medication", 3), rep("General PD Severity", 4), rep("Cognitive", 1) , rep("Disability", 9), rep("Sleep", 2), rep("Smell",1), rep("Gene expression",length(selected_genes)))
structure <- data.frame(node_names, community_structure)


dim(df_no_controls)

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
                        selected_genes)

structure <- subset(structure, structure$node_names %in% variables_to_scale)


# Loop through each variable and scale it
# for (variable in variables_to_scale) {
#   long_visit[[variable]] <- scale(long_visit[[variable]])
# }

#Check, patients are not filtered by idiopathic PD, all of them are included
table(long_visit$diagnosis_at_baseline) 


#3. Build multiple Linear Mixed Model (MariNET >= 2.0.0)

# Clinical variables are modeled on the rest of clinical variables and genes,
# and genes are modeled on clinical variables. Edge weight = mean of both t-values.
# Clinical-clinical and gene-gene interactions are set to 0 by the function.

clinical_variables <- setdiff(variables_to_scale, selected_genes)

symm_score_muted <- lmm_integration_analysis(
  as.data.frame(long_visit),               # Clinical data merged with gene expression
  variables_to_scale = clinical_variables, # Clinical variables (network nodes)
  molecular_features = selected_genes,     # Genes (network nodes)
  random_effects = "(1|sex) + (1|participant_id) + (1|age_at_baseline)"
)

# Same node order as community structure (clinical variables, then genes)
stopifnot(identical(colnames(symm_score_muted), structure$node_names))

spring_layout <- qgraph(symm_score_muted, layout = "spring", DoNotPlot = TRUE)$layout

# Plot new network (gene-clinical interactions only)

#png(here::here("results/transcriptomics/marinet/transcriptomics_clinical.png"), width = 2400, height = 1800, res = 300)  # 8x6 pulgadas a 300 dpi

qgraph(symm_score_muted, groups = structure$community_structure, layout=spring_layout, color=c("lightgreen", "lightblue", "#00BFFF","orange", "pink","grey"),
       labels = colnames(symm_score_muted)) 

#dev.off()


setwd(here::here("results/transcriptomics/marinet/"))

# Write output
score_df <- as.data.frame(symm_score_muted)
score_df<- cbind(row_names = rownames(score_df), score_df)
write.csv(score_df, "transcriptomic_clinical_score_matrix.csv")

# Paso 1: Calcular la suma de interacciones por nodo (suma por fila o columna, es igual porque es simétrico)
node_strength <- apply(abs(symm_score_muted[selected_genes, ]), 1, max)

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

el        <- make_edge_labels(submatrix, threshold = 2)
ec_labels <- make_edge_colors(el)
ec_lines  <- make_edge_colors_gradient(submatrix, exponent = 1.6)

fig_w <- 10     # inches
fig_h <- 7.5    # inches

plot_net <- function(mat, layout, groups = structure$community_structure) {
  qgraph(mat,
         vsize = 4,
         groups = groups,
         layout = layout,
         color = c("lightgreen", "lightblue", "#00BFFF", "orange", "pink", "grey"),
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

lay <- network_layout(submatrix)

png("transcriptomics_clinical_top10.png", width = fig_w, height = fig_h, units = "in", res = 300)
plot_net(submatrix, lay); dev.off()

svg("transcriptomics_clinical_top10.svg", width = fig_w, height = fig_h)
plot_net(submatrix, lay); dev.off()

pdf("transcriptomics_clinical_top10.pdf", width = fig_w, height = fig_h)
plot_net(submatrix, lay); dev.off()

score_output <- as.data.frame(submatrix)
score_output<- cbind(row_names = rownames(score_output), score_output)
write.csv(score_output, "transcriptomics_clinical_top10.csv")