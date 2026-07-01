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

metadata %>%
  dplyr::summarise(
    APOE_E4 = sum(has_known_APOE_E4_mutation_in_WGS == "Yes", na.rm = TRUE),
    LRRK2   = sum(has_known_LRRK2_mutation_in_WGS   == "Yes", na.rm = TRUE),
    GBA     = sum(has_known_GBA_mutation_in_WGS     == "Yes", na.rm = TRUE),
    SNCA    = sum(has_known_SNCA_mutation_in_WGS    == "Yes", na.rm = TRUE),
    Any_mutation = sum(has_known_PD_Mutation_in_WGS == "Yes", na.rm = TRUE),
    No_mutation  = sum(has_known_GBA_mutation_in_WGS     == "No" &
                         has_known_LRRK2_mutation_in_WGS   == "No" &
                         has_known_SNCA_mutation_in_WGS    == "No" &
                         has_known_APOE_E4_mutation_in_WGS == "No", na.rm = TRUE),
    Total = n()
  ) %>%
  pivot_longer(everything(), names_to = "Group", values_to = "N") %>%
  mutate(Pct = round(N / N[Group == "Total"] * 100, 1))


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
                "global_famhistory","PRS88","GBA_mutation", "LRRK2_mutation","SNCA_mutation",
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

multiple_visit %>%
  dplyr::summarise(
    APOE_E4 = sum(APOE_E4_mutation == "Yes", na.rm = TRUE),
    LRRK2   = sum(LRRK2_mutation   == "Yes", na.rm = TRUE),
    GBA     = sum(GBA_mutation     == "Yes", na.rm = TRUE),
    SNCA    = sum(SNCA_mutation    == "Yes", na.rm = TRUE),
    Any_mutation = sum(has_known_PD_Mutation_in_WGS == "Yes", na.rm = TRUE),
    No_mutation  = sum(GBA_mutation     == "No" &
                         LRRK2_mutation   == "No" &
                         SNCA_mutation    == "No" &
                         APOE_E4_mutation == "No", na.rm = TRUE),
    Total = n()
  ) %>%
  pivot_longer(everything(), names_to = "Group", values_to = "N") %>%
  mutate(Pct = round(N / N[Group == "Total"] * 100, 1))

length(unique(multiple_visit$participant_id))


# ---- Global PRS quartile groups (all participants) -------------------------
prs_global <- multiple_visit %>%
  distinct(participant_id, .keep_all = TRUE) %>%
  dplyr::select(participant_id, PRS88,
                APOE_E4_mutation, LRRK2_mutation, GBA_mutation, SNCA_mutation) %>%
  mutate(
    no_mutation = APOE_E4_mutation == "No" & LRRK2_mutation == "No" &
      GBA_mutation == "No"     & SNCA_mutation  == "No",
    prs_group = case_when(
      PRS88 <= q1_cutoff ~ "Q1 (Low PRS)",
      PRS88 >= q4_cutoff ~ "Q4 (High PRS)",
      TRUE               ~ "Q2-Q3"
    )
  ) %>%
  filter(prs_group != "Q2-Q3")  # keep only extremes

cat("Q1 (Low PRS):",  sum(prs_global$prs_group == "Q1 (Low PRS)"),  "participants\n")
cat("Q4 (High PRS):", sum(prs_global$prs_group == "Q4 (High PRS)"), "participants\n")

# ---- LMM on Q1 vs Q4 (all participants, no mutation stratification) --------
low_prs_ids  <- prs_global %>% filter(prs_group == "Q1 (Low PRS)")  %>% pull(participant_id)
high_prs_ids <- prs_global %>% filter(prs_group == "Q4 (High PRS)") %>% pull(participant_id)

low_prs_all  <- multiple_visit %>% filter(participant_id %in% low_prs_ids)
high_prs_all <- multiple_visit %>% filter(participant_id %in% high_prs_ids)

if(nrow(low_prs_all) >= 20){
  model_low_prs_all  <- lmm_analysis(low_prs_all,  variables_to_scale)
} else {
  model_low_prs_all  <- NULL
  cat("Skipping Low PRS — sample too small.\n")
}

if(nrow(high_prs_all) >= 20){
  model_high_prs_all <- lmm_analysis(high_prs_all, variables_to_scale)
} else {
  model_high_prs_all <- NULL
  cat("Skipping High PRS — sample too small.\n")
}

# ---- Plot: Q1 vs Q4 network comparison -------------------------------------
png("MariNET_PRS_Q1_vs_Q4.png",
    width = 2400, height = 1200, res = 300, type = "cairo", bg = "transparent")
layout(matrix(c(1, 2), nrow = 1))
par(mar = c(0, 0, 2, 0), cex.main = 1.2)

if(!is.null(model_low_prs_all)){
  edge_labels <- ifelse(abs(model_low_prs_all) > 3, round(model_low_prs_all, 2), NA)
  edge_label_colors <- ifelse(edge_labels > 0, "#33CC33", ifelse(edge_labels < 0, "red", "black"))
  qgraph(model_low_prs_all,
         title = "(A) Low PRS — Q1",
         layout = "spring", groups = structure_plot,
         color = c("lightgreen","lightblue","orange","pink","grey"),
         edge.labels = edge_labels, edge.label.cex = 1,
         edge.label.font = 2, edge.label.color = edge_label_colors,
         legend = FALSE, nodeNames = variables_to_scale,
         vsize = 7, label.cex = 1.2)
}

if(!is.null(model_high_prs_all)){
  edge_labels <- ifelse(abs(model_high_prs_all) > 3, round(model_high_prs_all, 2), NA)
  edge_label_colors <- ifelse(edge_labels > 0, "#33CC33", ifelse(edge_labels < 0, "red", "black"))
  qgraph(model_high_prs_all,
         title = "(B) High PRS — Q4",
         layout = "spring", groups = structure_plot,
         color = c("lightgreen","lightblue","orange","pink","grey"),
         edge.labels = edge_labels, edge.label.cex = 1,
         edge.label.font = 2, edge.label.color = edge_label_colors,
         legend = FALSE, nodeNames = variables_to_scale,
         vsize = 7, label.cex = 1.2)
}

layout(1)
dev.off()
cat("Saved: MariNET_PRS_Q1_vs_Q4.png\n")


# ---- Difference network: Q4 - Q1 -------------------------------------------
if(!is.null(model_low_prs_all) && !is.null(model_high_prs_all)){
  diff_prs_all <- model_high_prs_all - model_low_prs_all
  
  png("MariNET_PRS_difference.png",
      width = 1400, height = 1200, res = 300, type = "cairo", bg = "transparent")
  par(mar = c(0, 0, 2, 0))
  edge_labels <- ifelse(abs(diff_prs_all) > 3, round(diff_prs_all, 2), NA)
  edge_label_colors <- ifelse(edge_labels > 0, "#33CC33", ifelse(edge_labels < 0, "red", "black"))
  qgraph(diff_prs_all,
         title = "Difference network (Q4 - Q1)",
         layout = "spring", groups = structure_plot,
         color = c("lightgreen","lightblue","orange","pink","grey"),
         edge.labels = edge_labels, edge.label.cex = 1,
         edge.label.font = 2, edge.label.color = edge_label_colors,
         legend = FALSE, nodeNames = variables_to_scale,
         vsize = 7, label.cex = 1.2)
  dev.off()
  cat("Saved: MariNET_PRS_difference.png\n")
}
# ---- Difference network ----------------------------------------------------
# ---- Plot: Q1, Q4, Difference + Legend in one PNG -------------------------
png("MariNET_PRS_Q1_vs_Q4.png",
    width = 2400, height = 2400, res = 300, type = "cairo", bg = "transparent")

layout(matrix(c(1, 2,
                3, 4), nrow = 2, ncol = 2, byrow = TRUE))
par(mar = c(0, 0, 2, 0), cex.main = 1.2)

# Panel A: Low PRS (Q1)
if(!is.null(model_low_prs_all)){
  edge_labels       <- ifelse(abs(model_low_prs_all) > 3, round(model_low_prs_all, 2), NA)
  edge_label_colors <- ifelse(edge_labels > 0, "#33CC33", ifelse(edge_labels < 0, "red", "black"))
  qgraph(model_low_prs_all,
         title = "(A) Low PRS — Q4",
         layout = "spring", groups = structure_plot,
         color = c("lightgreen","lightblue","orange","pink","grey"),
         edge.labels = edge_labels, edge.label.cex = 1,
         edge.label.font = 2, edge.label.color = edge_label_colors,
         legend = FALSE, nodeNames = variables_to_scale,
         vsize = 7, label.cex = 1.2)
}

# Panel B: High PRS (Q4)
if(!is.null(model_high_prs_all)){
  edge_labels       <- ifelse(abs(model_high_prs_all) > 3, round(model_high_prs_all, 2), NA)
  edge_label_colors <- ifelse(edge_labels > 0, "#33CC33", ifelse(edge_labels < 0, "red", "black"))
  qgraph(model_high_prs_all,
         title = "(B) High PRS — Q1",
         layout = "spring", groups = structure_plot,
         color = c("lightgreen","lightblue","orange","pink","grey"),
         edge.labels = edge_labels, edge.label.cex = 1,
         edge.label.font = 2, edge.label.color = edge_label_colors,
         legend = FALSE, nodeNames = variables_to_scale,
         vsize = 7, label.cex = 1.2)
}

# Panel C: Difference (Q4 - Q1)
if(!is.null(model_low_prs_all) && !is.null(model_high_prs_all)){
  diff_prs_all      <- model_high_prs_all - model_low_prs_all
  edge_labels       <- ifelse(abs(diff_prs_all) > 2, round(diff_prs_all, 2), NA)
  edge_label_colors <- ifelse(edge_labels > 0, "#33CC33", ifelse(edge_labels < 0, "red", "black"))
  qgraph(diff_prs_all,
         title = "(C) Difference (Q1 - Q4)",
         layout = "spring", groups = structure_plot,
         color = c("lightgreen","lightblue","orange","pink","grey"),
         edge.labels = edge_labels, edge.label.cex = 1,
         edge.label.font = 2, edge.label.color = edge_label_colors,
         legend = FALSE, nodeNames = variables_to_scale,
         vsize = 7, label.cex = 1.2)
}

# Panel D: Legend (bottom right)
plot.new()
legend("center",
       legend = unique(structure_plot),
       fill   = c("lightgreen","lightblue","orange","pink","grey"),
       title  = "Community",
       cex    = 1.1,
       bty    = "n")

layout(1)
dev.off()
cat("Saved: MariNET_PRS_Q1_vs_Q4.png\n")



# ------------------------- PRS Stratification --------------------------------

# Compute quartile cutoffs per participant (one PRS value per person)
prs_per_participant <- multiple_visit %>%
  distinct(participant_id, .keep_all = TRUE) %>%
  dplyr::select(participant_id, PRS88)

q1_cutoff <- quantile(prs_per_participant$PRS88, 0.25, na.rm = TRUE)
q4_cutoff <- quantile(prs_per_participant$PRS88, 0.75, na.rm = TRUE)

cat("Q1 cutoff (<=):", q1_cutoff, "\n")
cat("Q4 cutoff (>=):", q4_cutoff, "\n")

# Get participant IDs for each extreme
low_prs_ids  <- prs_per_participant %>% filter(PRS88 <= q1_cutoff) %>% pull(participant_id)
high_prs_ids <- prs_per_participant %>% filter(PRS88 >= q4_cutoff) %>% pull(participant_id)

cat("Low PRS group (Q1):",  length(low_prs_ids),  "participants\n")
cat("High PRS group (Q4):", length(high_prs_ids), "participants\n")

# Filter all visits for those participants
low_prs_group  <- multiple_visit %>% filter(participant_id %in% low_prs_ids)
high_prs_group <- multiple_visit %>% filter(participant_id %in% high_prs_ids)

# ------------------------- LMM Analysis on PRS groups -----------------------

if(nrow(low_prs_group) >= 20){
  model_low_prs <- lmm_analysis(low_prs_group, variables_to_scale)
} else {
  model_low_prs <- NULL
  cat("Skipping Low PRS group due to low sample size.\n")
}

if(nrow(high_prs_group) >= 20){
  model_high_prs <- lmm_analysis(high_prs_group, variables_to_scale)
} else {
  model_high_prs <- NULL
  cat("Skipping High PRS group due to low sample size.\n")
}

df <- multiple_visit
df %>%
  distinct(participant_id, .keep_all = TRUE) %>%  # one row per participant
  dplyr::count(diagnosis_at_baseline)                      # count per diagnosis



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
# mut = "APOE_E4_mutation"
# mutated_group <- multiple_visit %>% filter(.data[[mut]] == "Yes")

# ------------------------- PRS Stratification per Mutation -------------------

# Compute quartile cutoffs once (per participant, not per row)
prs_per_participant <- multiple_visit %>%
  distinct(participant_id, .keep_all = TRUE) %>%
  dplyr::select(participant_id, PRS88)

q1_cutoff <- quantile(prs_per_participant$PRS88, 0.25, na.rm = TRUE)
q4_cutoff <- quantile(prs_per_participant$PRS88, 0.75, na.rm = TRUE)

cat("Q1 cutoff (<=):", q1_cutoff, "\n")
cat("Q4 cutoff (>=):", q4_cutoff, "\n")

# Store models per mutation x PRS stratum
models_prs <- list()

mutation_types <- c("APOE_E4_mutation", "LRRK2_mutation", "GBA_mutation", "SNCA_mutation")

for(mut in mutation_types){
  cat("\n--- Mutation:", mut, "---\n")
  
  mutated <- multiple_visit %>% filter(.data[[mut]] == "Yes")
  
  # Low PRS (Q1) within mutated group
  low_ids  <- mutated %>%
    distinct(participant_id, .keep_all = TRUE) %>%
    filter(PRS88 <= q1_cutoff) %>%
    pull(participant_id)
  
  high_ids <- mutated %>%
    distinct(participant_id, .keep_all = TRUE) %>%
    filter(PRS88 >= q4_cutoff) %>%
    pull(participant_id)
  
  low_group  <- mutated %>% filter(participant_id %in% low_ids)
  high_group <- mutated %>% filter(participant_id %in% high_ids)
  
  cat("Low PRS (Q1):",  length(low_ids),  "participants\n")
  cat("High PRS (Q4):", length(high_ids), "participants\n")
  
  if(nrow(low_group) >= 15){
    models_prs[[mut]][["low"]]  <- lmm_analysis(low_group,  variables_to_scale)
  } else {
    models_prs[[mut]][["low"]]  <- NULL
    cat("Skipping", mut, "Low PRS — sample too small.\n")
  }
  
  if(nrow(high_group) >= 15){
    models_prs[[mut]][["high"]] <- lmm_analysis(high_group, variables_to_scale)
  } else {
    models_prs[[mut]][["high"]] <- NULL
    cat("Skipping", mut, "High PRS — sample too small.\n")
  }
}

# ------------------------- Plot: one PNG per mutation ------------------------

for(mut in mutation_types){
  
  low_model  <- models_prs[[mut]][["low"]]
  high_model <- models_prs[[mut]][["high"]]
  
  # Skip if both are NULL
  if(is.null(low_model) && is.null(high_model)){
    cat("No models to plot for", mut, "\n")
    next
  }
  
  fname <- paste0("MariNET_PRS_", mut, ".png")
  png(fname, width = 2400, height = 1200, res = 300, type = "cairo", bg = "transparent")
  layout(matrix(c(1, 2), nrow = 1))
  par(mar = c(0, 0, 2, 0), cex.main = 1.2)
  
  # Plot Low PRS
  if(!is.null(low_model)){
    edge_labels <- ifelse(abs(low_model) > 3, round(low_model, 2), NA)
    edge_label_colors <- ifelse(edge_labels > 0, "#33CC33", ifelse(edge_labels < 0, "red", "black"))
    qgraph(low_model,
           title = paste0("(A) ", mut, " — Low PRS (Q1)"),
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
  } else {
    plot.new()
    title(paste0("(A) ", mut, " — Low PRS: insufficient n"))
  }
  
  # Plot High PRS
  if(!is.null(high_model)){
    edge_labels <- ifelse(abs(high_model) > 3, round(high_model, 2), NA)
    edge_label_colors <- ifelse(edge_labels > 0, "#33CC33", ifelse(edge_labels < 0, "red", "black"))
    qgraph(high_model,
           title = paste0("(B) ", mut, " — High PRS (Q4)"),
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
  } else {
    plot.new()
    title(paste0("(B) ", mut, " — High PRS: insufficient n"))
  }
  
  layout(1)
  dev.off()
  cat("Saved:", fname, "\n")
}

#########################################################################################

# ---- Rebuild long-format with "No mutation" as the reference group ----------
prs_mutations_clean <- multiple_visit %>%
  distinct(participant_id, .keep_all = TRUE) %>%
  dplyr::select(participant_id, PRS88,
                APOE_E4_mutation, LRRK2_mutation, GBA_mutation, SNCA_mutation) %>%
  mutate(no_mutation = case_when(
    APOE_E4_mutation == "No" & LRRK2_mutation == "No" &
      GBA_mutation == "No"     & SNCA_mutation  == "No" ~ TRUE,
    TRUE ~ FALSE
  ))

# Build one data frame per mutation: carriers vs no-mutation-at-all
prs_long_clean <- bind_rows(
  # APOE_E4
  prs_mutations_clean %>%
    filter(APOE_E4_mutation == "Yes" | no_mutation) %>%
    mutate(mutation = "APOE_E4",
           group    = ifelse(APOE_E4_mutation == "Yes", "Carrier", "No mutation")),
  # LRRK2
  prs_mutations_clean %>%
    filter(LRRK2_mutation == "Yes" | no_mutation) %>%
    mutate(mutation = "LRRK2",
           group    = ifelse(LRRK2_mutation == "Yes", "Carrier", "No mutation")),
  # GBA
  prs_mutations_clean %>%
    filter(GBA_mutation == "Yes" | no_mutation) %>%
    mutate(mutation = "GBA",
           group    = ifelse(GBA_mutation == "Yes", "Carrier", "No mutation")),
  # SNCA
  prs_mutations_clean %>%
    filter(SNCA_mutation == "Yes" | no_mutation) %>%
    mutate(mutation = "SNCA",
           group    = ifelse(SNCA_mutation == "Yes", "Carrier", "No mutation"))
)

# ---- Plot 1: Density — carrier vs no-mutation per mutation ------------------
p_density <- ggplot(prs_long_clean, aes(x = PRS88, fill = group, color = group)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  geom_vline(xintercept = q1_cutoff, linetype = "dashed", color = "red",       linewidth = 0.7) +
  geom_vline(xintercept = q4_cutoff, linetype = "dashed", color = "darkgreen", linewidth = 0.7) +
  facet_wrap(~ mutation, ncol = 2) +
  scale_fill_manual(values  = c("Carrier" = "#E04B4B", "No mutation" = "steelblue")) +
  scale_color_manual(values = c("Carrier" = "#E04B4B", "No mutation" = "steelblue")) +
  labs(title = "PRS88 Density: Mutation Carriers vs No Mutation",
       x = "PRS88", y = "Density", fill = NULL, color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

# ---- Plot 2: Boxplot — carrier vs no-mutation per mutation ------------------
p_boxplot <- ggplot(prs_long_clean, aes(x = group, y = PRS88, fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 1.2) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 0.7) +
  geom_hline(yintercept = q1_cutoff, linetype = "dashed", color = "red",       linewidth = 0.6) +
  geom_hline(yintercept = q4_cutoff, linetype = "dashed", color = "darkgreen", linewidth = 0.6) +
  facet_wrap(~ mutation, ncol = 2) +
  scale_fill_manual(values = c("Carrier" = "#E04B4B", "No mutation" = "steelblue")) +
  labs(title = "PRS88 Distribution: Mutation Carriers vs No Mutation",
       x = NULL, y = "PRS88") +
  theme_minimal() +
  theme(legend.position = "none")

# ---- Plot 3: Histogram — carrier vs no-mutation per mutation ----------------
p_hist <- ggplot(prs_long_clean, aes(x = PRS88, fill = group)) +
  geom_histogram(bins = 35, position = "identity", alpha = 0.55, color = "white") +
  geom_vline(xintercept = q1_cutoff, linetype = "dashed", color = "red",       linewidth = 0.7) +
  geom_vline(xintercept = q4_cutoff, linetype = "dashed", color = "darkgreen", linewidth = 0.7) +
  facet_wrap(~ mutation, ncol = 2) +
  scale_fill_manual(values = c("Carrier" = "#E04B4B", "No mutation" = "steelblue")) +
  labs(title = "PRS88 Histogram: Mutation Carriers vs No Mutation",
       x = "PRS88", y = "Count", fill = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

# ---- Plot 4: Sample sizes Q1/Q4 — carriers vs no-mutation ------------------
prs_quartile_clean <- prs_long_clean %>%
  mutate(prs_group = case_when(
    PRS88 <= q1_cutoff ~ "Q1 (Low)",
    PRS88 >= q4_cutoff ~ "Q4 (High)",
    TRUE               ~ "Q2-Q3 (Middle)"
  )) %>%
  dplyr::count(mutation, group, prs_group)

p_counts <- ggplot(prs_quartile_clean, aes(x = group, y = n, fill = prs_group)) +
  geom_col(position = "dodge", alpha = 0.85, color = "white") +
  geom_hline(yintercept = 20, linetype = "dashed", color = "black", linewidth = 0.7) +
  annotate("text", x = 0.5, y = 22, label = "n = 20", size = 3, hjust = 0) +
  facet_wrap(~ mutation, ncol = 2) +
  scale_fill_manual(values = c("Q1 (Low)"       = "#E04B4B",
                               "Q2-Q3 (Middle)" = "grey70",
                               "Q4 (High)"      = "#2E8B57")) +
  labs(title = "Sample Sizes by PRS Quartile: Carriers vs No Mutation",
       x = NULL, y = "N participants", fill = "PRS Group") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 15, hjust = 1))

# ---- Save ------------------------------------------------------------------
ggsave("PRS_density_carrier_vs_nomut.png",   p_density, width = 10, height = 7, dpi = 300)
ggsave("PRS_boxplot_carrier_vs_nomut.png",   p_boxplot, width = 10, height = 7, dpi = 300)
ggsave("PRS_histogram_carrier_vs_nomut.png", p_hist,    width = 10, height = 7, dpi = 300)
ggsave("PRS_counts_carrier_vs_nomut.png",    p_counts,  width = 10, height = 7, dpi = 300)

cat("All comparison plots saved.\n")








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
png("MariNET_PRS88_mutations.png", 
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
       color = c("lightgreen","lightblue","orange","magenta","pink","grey"),
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
       color = c("lightgreen","lightblue","orange","magenta","pink","grey"),
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
       color = c("lightgreen","lightblue","orange","magenta","pink","grey"),
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
       color = c("lightgreen","lightblue","orange","magenta","pink","grey"),
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



       