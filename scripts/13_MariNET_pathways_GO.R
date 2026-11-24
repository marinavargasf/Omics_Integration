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

setwd(here::here("results/pathways/")) #or results!

#1. Collect and filter data 
scores <- readRDS("pathwayScores_samplesFiltered_GO_collapse.rds")
rownames(scores) <- scores$term
scores <- scores[, !(colnames(scores) %in% c("ID", "term"))]

NES <- readRDS("nes_GO_collapse.RDS")

# traslation was previously relied upon as a loaded workspace object 
traslation <- readRDS("pathway_id_term_translation.rds")

NES <- NES%>%
  left_join(traslation, by = c("pathway" = "ID")) # Usa la columna ID como clave

NES$regulation <- ifelse(
  NES$NES > 0, "Positive regulated",
  ifelse(NES$NES < 0, "Negative regulated", NA)
)

NES$term <- gsub("[- ]", "_", NES$term)
#NES$pathway <- gsub("HALLMARK_", "", NES$pathway)


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
                               "SVM6T1" = 6,
                               "SVM6T1.1" = 6,
                               "SVM12T1" = 12,
                               "SVM12T1.1" = 12,
                               "SVM18T1" = 18,
                               "SVM18T1.1" = 18,
                               "SVM24T1" = 24,
                               "SVM24T1.1" = 24,
                               "SVM30T1" = 30,
                               "SVM30T1.1" = 30,
                               "SVM36T1" = 36,
                               "SVM36T1.1" = 36,
                               "SVM42T1" = 42,
                               "SVM42T1.1" = 42,
)

# Fail loudly rather than silently carrying NA visit_month into the merge/model steps
unmatched_time <- unique(selected_pathways$time[is.na(selected_pathways$visit_month)])
if (length(unmatched_time) > 0) {
  stop("Unmatched 'time' codes found in selected_pathways, these samples would get NA visit_month: ",
       paste(unmatched_time, collapse = ", "))
}

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

df_no_controls <- df_no_controls %>% dplyr::select( participant_id, visit_month, study, diagnosis_at_baseline, age_at_baseline, sex, race, case_control_other_at_baseline, global_famhistory,
                                            on_levodopa, on_dopamine_agonist, on_other_pd_medications, mds_updrs_part_i_summary_score, mds_updrs_part_ii_summary_score,
                                            mds_updrs_part_iii_summary_score, mds_updrs_part_iv_summary_score, moca_total_score, pdq39_mobility_score, pdq39_adl_score, pdq39_emotional_score, pdq39_stigma_score, pdq39_social_score, pdq39_cognition_score, pdq39_communication_score, pdq39_discomfort_score, 
                                            Schwad_ADL_score, REM_score, ess_sleepiness_score, upsit_total_score,
                                            all_of(pathway_names))

#pathway_names <- gsub("HALLMARK_", "", pathway_names)
pathway_names <- gsub("[- ]", "_", pathway_names)

# Select node names
node_names <- c("participant_id", "visit_month", "study", "diagnosis_at_baseline", "age_at_baseline", "sex", "race", "case_control_other_at_baseline", "global_famhistory"
                , "levodopa", "dopamine", "other_med", "UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4", "MOCA", "Mobility39", "ADL39","Emotional39", "Stigma39", "Social39", "Cognition39","Communication39","Discomfort39", "Schwad_ADL", "REM", "ESS", "UPSIT"
                ,pathway_names)
colnames(df_no_controls) <- node_names

# Induce community structure, and create dataset which contains this information
community_structure <- c(rep("Demographics", 9),rep("Medication", 3), rep("General PD Severity", 4), rep("Cognitive", 1) , rep("Disability", 9), rep("Sleep", 2), rep("Smell",1), rep("Gene expression",length(pathway_names)))
structure <- data.frame(node_names, community_structure)


###############################################################################
####################### Create stratified groups ##############################

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
                        pathway_names)

structure <- subset(structure, structure$node_names %in% variables_to_scale)

merged_df <- structure %>%
  left_join(NES, by = c("node_names" = "term"))  # Ajusta "pathway" al nombre correcto
merged_df$regulation[c(1:16)] <- merged_df$community_structure[c(1:16)]

# Loop through each variable and scale it
# for (variable in variables_to_scale) {
#   long_visit[[variable]] <- scale(long_visit[[variable]])
# }

#Check, patients are not filtered by idiopathic PD, all of them are included
table(long_visit$diagnosis_at_baseline) 


#3. Build multiple Linear Mixed Model (MariNET >= 2.0.0)

# Clinical variables are modeled on the rest of clinical variables and pathway scores,
# and pathway scores are modeled on clinical variables. Edge weight = mean of both t-values.
# Clinical-clinical and pathway-pathway interactions are set to 0 by the function.

clinical_variables <- setdiff(variables_to_scale, pathway_names)

symm_score_muted <- lmm_integration_analysis(
  as.data.frame(long_visit),               # Clinical data merged with pathway scores
  variables_to_scale = clinical_variables, # Clinical variables (network nodes)
  molecular_features = pathway_names,      # Pathway scores (network nodes)
  random_effects = "(1|sex) + (1|participant_id) + (1|age_at_baseline)"
)

# Same node order as community structure and regulation groups (clinical variables, then pathways)
stopifnot(identical(colnames(symm_score_muted), structure$node_names),
          identical(colnames(symm_score_muted), merged_df$node_names))

# Check against previous results (read only, nothing is overwritten)
old_file <- here::here("results/pathways/GO_clinical_score_matrix.csv")
if (file.exists(old_file)) {
  old_score <- read.csv(old_file, row.names = 1, check.names = FALSE)
  old_score <- as.matrix(old_score[, colnames(old_score) != "row_names"])
  rownames(old_score) <- colnames(old_score)
  stopifnot(setequal(colnames(old_score), colnames(symm_score_muted)))
  old_score <- old_score[colnames(symm_score_muted), colnames(symm_score_muted)]
  edges_old <- old_score[clinical_variables, pathway_names]
  edges_new <- symm_score_muted[clinical_variables, pathway_names]
  cat("Max |difference| vs previous matrix:", signif(max(abs(edges_new - edges_old)), 3), "\n")
  cat("Edges with |difference| > 0.01:", sum(abs(edges_new - edges_old) > 0.01), "of", length(edges_new), "\n")
}

# Layout now computed on the gene-clinical network (clinical-clinical edges are not estimated)
spring_layout <- qgraph(symm_score_muted, layout = "spring", DoNotPlot = TRUE)$layout

# ------------------------------------------------------------------------
# GO_GSEA plot (muted, clinical-clinical edges removed) -- PNG, SVG, PDF
# ------------------------------------------------------------------------

png(here::here("results/pathways/GO_GSEA_MariNET2.png"), width = 2400, height = 1800, res = 300)

qgraph(
  symm_score_muted,
  groups = merged_df$regulation,
  layout = spring_layout,
  color = c("lightgreen", "lightblue", "orange", "coral", "#00BFFF", "pink", "grey"),
  labels = colnames(symm_score_muted),
  label.cex = 0.3,
  label.scale = FALSE,
  label.scale.equal = TRUE,
  vsize = 2,
  posCol = "#B0E57C", negCol =  "#FF7F7F"# first color = negative edges, second = positive edges
)
dev.off()

svg(here::here("results/pathways/GO_GSEA_MariNET2.svg"), width = 2400/300, height = 1800/300)

qgraph(
  symm_score_muted,
  groups = merged_df$regulation,
  layout = spring_layout,
  color = c("lightgreen", "lightblue", "orange", "coral", "#00BFFF", "pink", "grey"),
  labels = colnames(symm_score_muted),
  label.cex = 0.3,
  label.scale = FALSE,
  label.scale.equal = TRUE,
  vsize = 2,
  posCol = "#B0E57C", negCol =  "#FF7F7F"
)
dev.off()

pdf(here::here("results/pathways/GO_GSEA_MariNET2.pdf"), width = 2400/300, height = 1800/300)

qgraph(
  symm_score_muted,
  groups = merged_df$regulation,
  layout = spring_layout,
  color = c("lightgreen", "lightblue", "orange", "coral", "#00BFFF", "pink", "grey"),
  labels = colnames(symm_score_muted),
  label.cex = 0.3,
  label.scale = FALSE,
  label.scale.equal = TRUE,
  vsize = 2,
  posCol = "#B0E57C", negCol =  "#FF7F7F"
)
dev.off()

symm_score_muted_tograph <- symm_score_muted

symm_score_muted_tograph[abs(symm_score_muted_tograph)<0.8] <-0

# ------------------------------------------------------------------------
# GO_GSEA_clean plot (muted + thresholded |t|>=1) -- PNG, SVG, PDF
# ------------------------------------------------------------------------

png(here::here("results/pathways/GO_GSEA_clean_MariNET2.png"), width = 2400, height = 1800, res = 300)

qgraph(
  symm_score_muted_tograph,
  groups = merged_df$regulation,
  layout = spring_layout,
  color = c("lightgreen", "lightblue", "orange", "coral", "#00BFFF", "pink", "grey"),
  labels = colnames(symm_score_muted_tograph),
  label.cex = 0.3,
  label.scale = FALSE,
  label.scale.equal = TRUE,
  vsize = 2,
  posCol = "#B0E57C", negCol =  "#FF7F7F"# first color = negative edges, second = positive edges
)

dev.off()

svg(here::here("results/pathways/GO_GSEA_clean_MariNET2.svg"), width = 2400/300, height = 1800/300)

qgraph(
  symm_score_muted_tograph,
  groups = merged_df$regulation,
  layout = spring_layout,
  color = c("lightgreen", "lightblue", "orange", "coral", "#00BFFF", "pink", "grey"),
  labels = colnames(symm_score_muted_tograph),
  label.cex = 0.3,
  label.scale = FALSE,
  label.scale.equal = TRUE,
  vsize = 2,
  posCol = "#B0E57C", negCol =  "#FF7F7F"
)

dev.off()

pdf(here::here("results/pathways/GO_GSEA_clean_MariNET2.pdf"), width = 2400/300, height = 1800/300)

qgraph(
  symm_score_muted_tograph,
  groups = merged_df$regulation,
  layout = spring_layout,
  color = c("lightgreen", "lightblue", "orange", "coral", "#00BFFF", "pink", "grey"),
  labels = colnames(symm_score_muted_tograph),
  label.cex = 0.3,
  label.scale = FALSE,
  label.scale.equal = TRUE,
  vsize = 2,
  posCol = "#B0E57C", negCol =  "#FF7F7F"
)

dev.off()


setwd(here::here("results/pathways/"))

# Write output
score_df <- as.data.frame(symm_score_muted)
score_df<- cbind(row_names = rownames(score_df), score_df)
write.csv(score_df, "GO_clinical_score_matrix_MariNET2.csv")

# # Paso 1: Calcular la suma de interacciones por nodo (suma por fila o columna, es igual porque es simétrico)
# node_strength <- rowSums(abs(symm_score_muted[pathway_names,]))
# 
# # Paso 2: Convertirlo a dataframe
# strength_df <- data.frame(
#   node = names(node_strength),
#   strength = node_strength
# )
# 
# # Paso 4: Filtrar solo "protein expression" y seleccionar los 10 más conectados
# top10_genes <- rownames(strength_df[order(-strength_df$strength), ][1:10, ])
# 
# 
# # Resultado: nodos de proteína más conectados
# print(top10_genes)
# 
# selected_vars <- c("UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4", "MOCA", "Mobility39", 
#                    "ADL39", "Emotional39", "Stigma39", "Social39", "Cognition39",
#                    "Communication39", "Discomfort39", "Schwad_ADL", "ESS", "UPSIT",
#                    top10_genes)
# 
# submatrix <- symm_score_muted[selected_vars, selected_vars]
# structure <- subset(structure, structure$node_names %in% selected_vars)
# 
# merged_df <- structure %>%
#   left_join(NES, by = c("node_names" = "pathway"))  # Ajusta "pathway" al nombre correcto
# 
# # merged_df$regulation[c(1:16)] <- merged_df$community_structure[c(1:16)]
# 
# # Open a PNG graphics device
# png("pathways_clinical_top10.png", width = 1800, height = 1000, res = 500)
# 
# # Generate the qgraph plot
# qgraph(
#   submatrix,
#   groups = merged_df$community_structure,
#   layout = spring_layout,
#   color = c("lightgreen", "lightblue","coral", "orange",  "#00BFFF", "pink", "grey"),
#   legend.cex = 0.2,
#   labels = colnames(submatrix),
#   label.cex = 0.3,
#   label.scale = FALSE,
#   label.scale.equal = TRUE,
#   vsize = 2,
#   posCol = "#B0E57C", negCol =  "#FF7F7F"# first color = negative edges, second = positive edges
# )
# 
#
# # Close the graphics device
# dev.off()
# 
# score_output <- as.data.frame(submatrix)
# score_output<- cbind(row_names = rownames(score_output), score_output)
# write.csv(score_output, "pathways_clinical_top10.csv")