library(MariNET)
library(qgraph)
library(dplyr)

# ------------------------- Data Loading ---------------------------

source("~/Documents/Parkinson/clinical_metadata_v4_2023/DataReading_v4.R")

source("~/Documents/Parkinson/data_processing.R") 
setwd("~/Documents/Omics_Integration/Genomics_analysis")
polygenic <- read.csv("prs.csv")
colnames(polygenic)[1] <- "participant_id"
metadata <- merge(polygenic, metadata, by = "participant_id", all.x = TRUE)

# Merge clinical_data and metadata
master_matrix <- merge(clinical_data, metadata, by = "participant_id", all.x = TRUE)

# Select relevant columns
master_matrix <- master_matrix %>% select(
  participant_id, visit_month, study, diagnosis_at_baseline, age_at_baseline, sex, race, 
  case_control_other_at_baseline, global_famhistory, PRS88,
  has_known_GBA_mutation_in_WGS, has_known_LRRK2_mutation_in_WGS,
  has_known_SNCA_mutation_in_WGS, has_known_APOE_E4_mutation_in_WGS,
  has_known_PD_Mutation_in_WGS,
  mds_updrs_part_i_summary_score, mds_updrs_part_ii_summary_score, 
  mds_updrs_part_iii_summary_score, mds_updrs_part_iv_summary_score, 
  code_upd2hy_hoehn_and_yahr_stage, moca_total_score,
  pdq39_mobility_score, pdq39_adl_score, pdq39_emotional_score, pdq39_stigma_score, 
  pdq39_social_score, pdq39_cognition_score, pdq39_communication_score,
  pdq39_discomfort_score, Schwad_ADL_score, REM_score, ess_sleepiness_score, 
  upsit_total_score
)

# Keep only PD patients
master_matrix <- master_matrix %>% filter(diagnosis_at_baseline != "No PD Nor Other Neurological Disorder")

# Count unique participants
length(unique(master_matrix$participant_id))



# Rename columns
node_names <- c("participant_id", "visit_month", "study", "diagnosis_at_baseline", 
                "age_at_baseline", "sex", "race", "case_control_other_at_baseline", 
                "global_famhistory", "PRS88","GBA_mutation", "LRRK2_mutation","SNCA_mutation",
                "APOE_E4_mutation", "has_known_PD_Mutation_in_WGS",
                "UPDRS1", "UPDRS2","UPDRS3", "UPDRS4", "HOEHN_Stage", "MOCA", 
                "Mobility39", "ADL39","Emotional39", "Stigma39", "Social39", 
                "Cognition39", "Communication39","Discomfort39", "Schwad_ADL", 
                "REM", "ESS", "UPSIT")

colnames(master_matrix) <- node_names

# Define community structure
community_structure <- c(
  rep("Demographics", 9),
  rep("Polygenic Risk Score",1),
  rep("Genomics", 5),
  rep("General PD Severity", 5),
  rep("Cognitive", 1),
  rep("Disability", 9),
  rep("Sleep", 2),
  rep("Smell", 1)
)
structure <- data.frame(node_names, community_structure)

# Filter rows with <=7 missing values
multiple_visit <- master_matrix %>% dplyr::filter(rowSums(is.na(.)) <= 7)

# Filter columns with <70% missing
umbral <- 0.7
missing_pct <- colMeans(is.na(multiple_visit))
columnas_filtradas <- names(missing_pct[missing_pct < umbral])
multiple_visit <- multiple_visit[, columnas_filtradas]


length(unique(multiple_visit$participant_id))

df <- multiple_visit
df %>%
  distinct(participant_id, .keep_all = TRUE) %>%  # one row per participant
  count(diagnosis_at_baseline)                      # count per diagnosis


# Update structure to keep only filtered columns
structure <- subset(structure, node_names %in% columnas_filtradas)

# Define variables to scale
variables_to_scale <- c("UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4", "MOCA", "Mobility39",
                        "ADL39", "Emotional39", "Stigma39", "Social39", "Cognition39",
                        "Communication39", "Discomfort39", "Schwad_ADL", "ESS", "UPSIT")

# Define structure_plot for qgraph
structure_plot <- structure %>%
  filter(node_names %in% variables_to_scale) %>%
  pull(community_structure)

# Check for mutations
mutation_vars <- c("has_known_GBA_mutation_in_WGS", "has_known_LRRK2_mutation_in_WGS",
                   "has_known_SNCA_mutation_in_WGS", "has_known_APOE_E4_mutation_in_WGS",
                   "has_known_PD_Mutation_in_WGS")

mutation_summary <- lapply(mutation_vars, function(var){
  table(metadata[[var]], useNA = "ifany")
})

names(mutation_summary) <- mutation_vars
mutation_summary

# Patients with NO known mutation
no_mutation_group <- multiple_visit %>%
  filter(GBA_mutation == "No" &
           LRRK2_mutation == "No" &
           SNCA_mutation == "No" &
           APOE_E4_mutation == "No")

#LMM analysis on no mutated 
if(nrow(no_mutation_group) >= 20){
  model_no_mut <- lmm_analysis(no_mutation_group, variables_to_scale)
} else {
  model_no_mut <- NULL
  cat("Skipping No Mutation group due to low sample size.\n")
}
# ------------------------- Run MariNET for APOE and LRRK2 -------------------
mutation_types <- c("APOE_E4_mutation", "LRRK2_mutation", "GBA_mutation", "SNCA_mutation")
#mutation_types <- c("GBA_mutation")


models <- list()

for(mut in mutation_types){
  cat("Running MariNET for mutation:", mut, "\n")
  
  mutated_group <- multiple_visit %>% filter(.data[[mut]] == "Yes")
  
  if(nrow(mutated_group) < 20){
    cat("Skipping", mut, "due to low sample size.\n")
    next
  }
  
  model_mut <- lmm_analysis(mutated_group, variables_to_scale)
  models[[mut]] <- model_mut
}

# Compute difference network
#difference_mut <- models[["LRRK2_mutation"]] - models[["GBA_mutation"]]

# ------------------------- Plot Networks with original style -----------------
png("MariNET_mutations.png", 
    width = 2400, height = 2000, res = 300, type = "cairo", bg = "transparent")

# 2x2 layout: Top-left (APOE), Top-right (LRRK2), Bottom-left (Difference), Bottom-right (Legend)
layout_matrix <- matrix(c(1,2,3,4), nrow=2, ncol=2, byrow=TRUE)
layout(layout_matrix)
par(mar=c(0,0,2,0), cex.main=1.2)  # Margins and title font

# ------------------------- Plot 1: APOE_E4 MariNET ----------------------------
edge_labels <- ifelse(abs(models[["APOE_E4_mutation"]]) > 3, round(models[["APOE_E4_mutation"]], 2), NA)
# Optional: print specific interaction like before
edge_labels[5,14] = edge_labels[14,5] <- round(models[["APOE_E4_mutation"]][5,14], 3)
edge_label_colors <- ifelse(edge_labels > 0, "#33CC33", ifelse(edge_labels < 0, "red", "black"))

qgraph(models[["APOE_E4_mutation"]],
       title = "(A) APOE_E4 MariNET",
       layout = "spring",
       groups = structure_plot,
       color = c("lightgreen","lightblue","orange","pink","grey"),
       edge.labels = edge_labels,
       edge.label.cex = 1,
       edge.label.font = 2,
       edge.label.color = edge_label_colors,
       legend = FALSE,
       nodeNames = variables_to_scale,
       vsize = 7,
       label.cex = 1.2)

# ------------------------- Plot 2: LRRK2 MariNET ----------------------------
edge_labels <- ifelse(abs(models[["LRRK2_mutation"]]) > 3, round(models[["LRRK2_mutation"]], 2), NA)
edge_label_colors <- ifelse(edge_labels > 0, "#33CC33", ifelse(edge_labels < 0, "red", "black"))

qgraph(models[["LRRK2_mutation"]],
       title = "(B) LRRK2 MariNET",
       layout = "spring",
       groups = structure_plot,
       color = c("lightgreen","lightblue","orange","pink","grey"),
       edge.labels = edge_labels,
       edge.label.cex = 1,
       edge.label.font = 2,
       edge.label.color = edge_label_colors,
       legend = FALSE,
       nodeNames = variables_to_scale,
       vsize = 7,
       label.cex = 1.2)

# ------------------------- Plot 3: Difference network ------------------------
# edge_labels <- ifelse(abs(difference_mut) > 5, round(difference_mut, 2), NA)
# edge_label_colors <- ifelse(edge_labels > 0, "#33CC33", ifelse(edge_labels < 0, "red", "black"))
# 
# qgraph(difference_mut,
#        title = "(C) Difference (APOE_E4 - LRRK2)",
#        layout = "spring",
#        groups = structure_plot,
#        color = c("lightgreen","lightblue","orange","pink","grey"),
#        edge.labels = edge_labels,
#        edge.label.cex = 1,
#        edge.label.font = 2,
#        edge.label.color = edge_label_colors,
#        legend = FALSE,
#        nodeNames = variables_to_scale,
#        vsize = 7,
#        label.cex = 1.2)
edge_labels <- ifelse(abs(models[["GBA_mutation"]]) > 3, round(models[["GBA_mutation"]], 2), NA)
edge_label_colors <- ifelse(edge_labels > 0, "#33CC33", ifelse(edge_labels < 0, "red", "black"))

qgraph(models[["GBA_mutation"]],
       title = "(C) GBA MariNET",
       layout = "spring",
       groups = structure_plot,
       color = c("lightgreen","lightblue","orange","pink","grey"),
       edge.labels = edge_labels,
       edge.label.cex = 1,
       edge.label.font = 2,
       edge.label.color = edge_label_colors,
       legend = FALSE,
       nodeNames = variables_to_scale,
       vsize = 7,
       label.cex = 1.2)

# -----------------
# ------------------------- Plot 4: Legend ------------------------------------
# plot(1, type="n", axes=FALSE, xlab="", ylab="")
# legend("center",
#        legend = unique(structure_plot),
#        fill = c("lightgreen","lightblue","orange","pink","grey"),
#        cex = 1.1,
#        bty = "n")
edge_labels <- ifelse(abs(model_no_mut) > 3, round(model_no_mut, 2), NA)
edge_label_colors <- ifelse(edge_labels > 0, "#33CC33", ifelse(edge_labels < 0, "red", "black"))

qgraph(model_no_mut,
       title = "(C) No mutation",
       layout = "spring",
       groups = structure_plot,
       color = c("lightgreen","lightblue","orange","pink","grey"),
       edge.labels = edge_labels,
       edge.label.cex = 1,
       edge.label.font = 2,
       edge.label.color = edge_label_colors,
       legend = FALSE,
       nodeNames = variables_to_scale,
       vsize = 7,
       label.cex = 1.2)


# Reset layout
layout(1)
dev.off()



#PLOT difference between mutations and non mutated
# ---- Helper: safe edge labels and colors -----------------------------------
make_edge_labels <- function(model, threshold = 3){
  el <- matrix(NA, nrow = nrow(model), ncol = ncol(model))
  el[abs(model) > threshold] <- round(model[abs(model) > threshold], 2)
  el
}

make_edge_colors <- function(el){
  ec <- matrix("black", nrow = nrow(el), ncol = ncol(el))
  ec[!is.na(el) & el > 0] <- "#33CC33"
  ec[!is.na(el) & el < 0] <- "red"
  ec
}

# ---- Plot loop -------------------------------------------------------------
for(mut in mutation_types){
  
  model_mut <- models[[mut]]
  if(is.null(model_mut)){
    cat("No model for", mut, "— skipping plot.\n")
    next
  }
  
  diff_model <- NULL
  if(!is.null(model_no_mut)){
    diff_model <- model_mut - model_no_mut
  }
  
  mut_label <- gsub("_mutation", "", mut)
  fname     <- paste0("MariNET_", mut, ".png")
  
  png(fname, width = 2400, height = 2400, res = 300, type = "cairo", bg = "transparent")
  layout(matrix(c(1, 2,
                  3, 4), nrow = 2, ncol = 2, byrow = TRUE))
  par(mar = c(0, 0, 2, 0), cex.main = 1.2)
  
  # Panel A: Mutation network
  el <- make_edge_labels(model_mut)
  ec <- make_edge_colors(el)
  qgraph(model_mut,
         title = paste0("(A) ", mut_label, " carriers"),
         layout = "spring", groups = structure_plot,
         color = c("lightgreen","lightblue","orange","pink","grey"),
         edge.labels = el, edge.label.cex = 1,
         edge.label.font = 2, edge.label.color = ec,
         legend = FALSE, nodeNames = variables_to_scale,
         vsize = 7, label.cex = 1.2)
  
  # Panel B: No mutation network
  if(!is.null(model_no_mut)){
    el <- make_edge_labels(model_no_mut)
    ec <- make_edge_colors(el)
    qgraph(model_no_mut,
           title = "(B) No mutation",
           layout = "spring", groups = structure_plot,
           color = c("lightgreen","lightblue","orange","pink","grey"),
           edge.labels = el, edge.label.cex = 1,
           edge.label.font = 2, edge.label.color = ec,
           legend = FALSE, nodeNames = variables_to_scale,
           vsize = 7, label.cex = 1.2)
  } else {
    plot.new(); title("(B) No mutation — unavailable")
  }
  
  # Panel C: Difference
  if(!is.null(diff_model)){
    el <- make_edge_labels(diff_model,threshold = 3)
    ec <- make_edge_colors(el)
    qgraph(diff_model,
           title = paste0("(C) Difference (", mut_label, " - No mutation)"),
           layout = "spring", groups = structure_plot,
           color = c("lightgreen","lightblue","orange","pink","grey"),
           edge.labels = el, edge.label.cex = 1,
           edge.label.font = 2, edge.label.color = ec,
           legend = FALSE, nodeNames = variables_to_scale,
           vsize = 7, label.cex = 1.2)
  } else {
    plot.new(); title("(C) Difference — unavailable")
  }
  
  # Panel D: Legend
  plot.new()
  legend("center",
         legend = unique(structure_plot),
         fill   = c("lightgreen","lightblue","orange","pink","grey"),
         title  = "Community",
         cex    = 1.1,
         bty    = "n")
  
  layout(1)
  dev.off()
  cat("Saved:", fname, "\n")
}
# As a summary dataframe
sample_size_df <- data.frame(
  group = c("No mutation", mutation_types),
  n_visits = c(
    nrow(no_mutation_group),
    sapply(mutation_types, function(mut) nrow(multiple_visit %>% filter(.data[[mut]] == "Yes")))
  ),
  n_participants = c(
    length(unique(no_mutation_group$participant_id)),
    sapply(mutation_types, function(mut) length(unique(multiple_visit %>% filter(.data[[mut]] == "Yes") %>% pull(participant_id))))
  )
)

print(sample_size_df)


# ------------------------- Subsampled analysis for comparability -------------

set.seed(111)  # for reproducibility

models_subsampled <- list()

for(mut in mutation_types){
  mutated_group <- multiple_visit %>% filter(.data[[mut]] == "Yes")
  n_mut <- length(unique(mutated_group$participant_id))
  
  if(n_mut < 20){
    cat("Skipping", mut, "due to low sample size.\n")
    next
  }
  
  # Subsample no-mutation group to match mutation group size
  no_mut_participants <- unique(no_mutation_group$participant_id)
  sampled_participants <- sample(no_mut_participants, n_mut, replace = FALSE)
  no_mutation_subsampled <- no_mutation_group %>% 
    filter(participant_id %in% sampled_participants)
  
  cat(mut, "— mutation n:", n_mut, 
      "| no-mutation subsampled n:", 
      length(unique(no_mutation_subsampled$participant_id)), "\n")
  
  # Fit models on matched sample sizes
  model_mut        <- lmm_analysis(mutated_group, variables_to_scale)
  model_no_mut_sub <- lmm_analysis(no_mutation_subsampled, variables_to_scale)
  
  models_subsampled[[mut]] <- list(
    mutation   = model_mut,
    no_mutation = model_no_mut_sub,
    difference  = model_mut - model_no_mut_sub
  )
}

# ------------------------- Sample sizes per subgroup -------------------------
sample_size_df <- data.frame(
  group = c("No mutation", mutation_types),
  n_visits = c(
    nrow(no_mutation_group),
    sapply(mutation_types, function(mut) nrow(multiple_visit %>% filter(.data[[mut]] == "Yes")))
  ),
  n_participants = c(
    length(unique(no_mutation_group$participant_id)),
    sapply(mutation_types, function(mut) length(unique(multiple_visit %>% filter(.data[[mut]] == "Yes") %>% pull(participant_id))))
  )
)
print(sample_size_df)

# ---- Plot loop with subsampled no-mutation group ---------------------------
set.seed(111)

for(mut in mutation_types){
  
  model_mut <- models[[mut]]
  if(is.null(model_mut)){
    cat("No model for", mut, "— skipping plot.\n")
    next
  }
  
  # Subsample no-mutation group to match mutation group participant size
  n_mut <- sample_size_df$n_participants[sample_size_df$group == mut]
  no_mut_participants <- unique(no_mutation_group$participant_id)
  
  if(n_mut < 20){
    cat("Skipping", mut, "— mutation group too small for subsampling.\n")
    next
  }
  
  sampled_participants <- sample(no_mut_participants, n_mut, replace = FALSE)
  no_mutation_subsampled <- no_mutation_group %>%
    filter(participant_id %in% sampled_participants)
  
  cat(mut, "— mutation n:", n_mut,
      "| no-mutation subsampled n:", 
      length(unique(no_mutation_subsampled$participant_id)), "\n")
  
  model_no_mut_sub <- lmm_analysis(no_mutation_subsampled, variables_to_scale)
  
  diff_model <- NULL
  if(!is.null(model_no_mut_sub)){
    diff_model <- model_mut - model_no_mut_sub
  }
  
  mut_label <- gsub("_mutation", "", mut)
  fname     <- paste0("MariNET_", mut, "_subsampled.png")
  
  png(fname, width = 2400, height = 2400, res = 300, type = "cairo", bg = "transparent")
  layout(matrix(c(1, 2,
                  3, 4), nrow = 2, ncol = 2, byrow = TRUE))
  par(mar = c(0, 0, 2, 0), cex.main = 1.2)
  
  # Panel A: Mutation network
  el <- make_edge_labels(model_mut)
  ec <- make_edge_colors(el)
  qgraph(model_mut,
         title = paste0("(A) ", mut_label, " carriers"),
         layout = "spring", groups = structure_plot,
         color = c("lightgreen","lightblue","orange","pink","grey"),
         edge.labels = el, edge.label.cex = 1,
         edge.label.font = 2, edge.label.color = ec,
         legend = FALSE, nodeNames = variables_to_scale,
         vsize = 7, label.cex = 1.2)
  
  # Panel B: Subsampled no-mutation network
  el <- make_edge_labels(model_no_mut_sub)
  ec <- make_edge_colors(el)
  qgraph(model_no_mut_sub,
         title = paste0("(B) No mutation (subsampled n=", n_mut, ")"),
         layout = "spring", groups = structure_plot,
         color = c("lightgreen","lightblue","orange","pink","grey"),
         edge.labels = el, edge.label.cex = 1,
         edge.label.font = 2, edge.label.color = ec,
         legend = FALSE, nodeNames = variables_to_scale,
         vsize = 7, label.cex = 1.2)
  
  # Panel C: Difference
  if(!is.null(diff_model)){
    el <- make_edge_labels(diff_model, threshold = 3)
    ec <- make_edge_colors(el)
    qgraph(diff_model,
           title = paste0("(C) Difference (", mut_label, " - No mutation)"),
           layout = "spring", groups = structure_plot,
           color = c("lightgreen","lightblue","orange","pink","grey"),
           edge.labels = el, edge.label.cex = 1,
           edge.label.font = 2, edge.label.color = ec,
           legend = FALSE, nodeNames = variables_to_scale,
           vsize = 7, label.cex = 1.2)
  } else {
    plot.new(); title("(C) Difference — unavailable")
  }
  
  # Panel D: Legend
  plot.new()
  legend("center",
         legend = unique(structure_plot),
         fill   = c("lightgreen","lightblue","orange","pink","grey"),
         title  = "Community",
         cex    = 1.1,
         bty    = "n")
  
  layout(1)
  dev.off()
  cat("Saved:", fname, "\n")
}
       
