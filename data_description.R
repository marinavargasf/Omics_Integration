library(dplyr)

# ---- Transcriptomics -----------------------------------------------------
metadata_tx <- readRDS("~/Documents/Omics_Integration/Transcriptomic_analysis/DEA2026/data/confoundingDataFrame2_full.rds")
metadata_tx$participant_id <- sub("-[^-]+$", "", rownames(metadata_tx))
participants_tx <- metadata_tx %>%
  group_by(participant_id) %>%
  slice(1) %>%
  ungroup()

cat("Transcriptomic cohort — N participants:", nrow(participants_tx), "\n")
cat("Mean age:", round(mean(participants_tx$age_at_baseline, na.rm = TRUE), 1),
    "± SD", round(sd(participants_tx$age_at_baseline, na.rm = TRUE), 1), "\n")
cat("Sex distribution:\n"); print(table(participants_tx$sex))
cat("Race distribution:\n"); print(table(participants_tx$race))
cat("Case/control:\n"); print(table(participants_tx$case_control))

case_control_by_time <- metadata_tx %>%
  distinct(participant_id, time, visit_month, case_control) %>%
  count(visit_month, time, case_control) %>%
  tidyr::pivot_wider(names_from = case_control, values_from = n, values_fill = 0) %>%
  rename(Cases = case, Controls = control) %>%
  mutate(Total = Cases + Controls) %>%
  arrange(visit_month)
print(case_control_by_time)

participants_tx_ids <- unique(participants_tx$participant_id)
length(participants_tx_ids)

# ---- Shared proteomics metadata -------------------------------------------
metadata_prot <- readRDS(file = "~/Documents/Parkinson/proteomic/metadata.rds")

# ---- CSF -------------------------------------------------------------------
csf_dir <- "~/Documents/Parkinson/proteomic/Targeted olink/proteomics-CSF-PPEA-D03/"
prefix_csf <- "releases_2023_v4release_1027_proteomics-CSF-PPEA-D03_"

samples_CSF_raw <- read.csv(paste0(csf_dir, prefix_csf, "proteomics_CSF-PPEA-D03_samples.csv"))
cat("CSF — original samples (raw sample sheet):", nrow(samples_CSF_raw), "\n")

panel_dir_csf <- paste0(csf_dir, "protein-expression/", prefix_csf, "olink-explore_protein-expression_CSF-PPEA-D03_matrix_")
cardiometa_csf   <- read.csv(paste0(panel_dir_csf, "cardiometabolic.csv"))
inflammatory_csf <- read.csv(paste0(panel_dir_csf, "inflammation.csv"))
neurology_csf    <- read.csv(paste0(panel_dir_csf, "neurology.csv"))
oncology_csf     <- read.csv(paste0(panel_dir_csf, "oncology.csv"))

proteomic_raw_csf <- rbind(cardiometa_csf, inflammatory_csf, neurology_csf, oncology_csf)
rownames(proteomic_raw_csf) <- proteomic_raw_csf[[1]]
proteomic_raw_csf <- proteomic_raw_csf[, -1]
cat("CSF — samples in protein-expression matrix:", ncol(proteomic_raw_csf), "\n")

samples_CSF_full <- merge(samples_CSF_raw, metadata_prot[, c(1, 11, 13, 14)])
colnames(samples_CSF_full)[colnames(samples_CSF_full) == "case_control_other_at_baseline"] <- "condition"
samples_CSF <- samples_CSF_full %>% filter(condition %in% c("Case", "Control"))



participants_csf <- samples_CSF %>% group_by(participant_id) %>% slice(1) %>% ungroup()
participants_csf_ids <- unique(participants_csf$participant_id)

cat("CSF — N participants:", nrow(participants_csf), "\n")
cat("Mean age:", round(mean(participants_csf$age_at_baseline, na.rm = TRUE), 1),
    "± SD", round(sd(participants_csf$age_at_baseline, na.rm = TRUE), 1), "\n")
cat("Sex distribution:\n"); print(table(participants_csf$sex))
cat("Condition distribution:\n"); print(table(participants_csf$condition))

case_control_csf_by_time <- samples_CSF %>%
  count(visit_month, condition) %>%
  tidyr::pivot_wider(names_from = condition, values_from = n, values_fill = 0) %>%
  rename(Cases = Case, Controls = Control) %>%
  mutate(Total = Cases + Controls) %>%
  arrange(visit_month)
print(case_control_csf_by_time)

cat("\nCSF — total samples:", sum(case_control_csf_by_time$Total), "\n")
cat("CSF — total cases:", sum(case_control_csf_by_time$Cases), "\n")
cat("CSF — total controls:", sum(case_control_csf_by_time$Controls), "\n")

# ---- Plasma ------------------------------------------------------------------
pla_dir <- "~/Documents/Parkinson/proteomic/Targeted olink/proteomics-PLA-PPEA-D03/"
prefix_pla <- "releases_2023_v4release_1027_proteomics-PLA-PPEA-D03_"

samples_plasma_raw <- read.csv(paste0(pla_dir, prefix_pla, "proteomics_PLA-PPEA-D03_samples.csv"))
cat("PLA — original samples (raw sample sheet):", nrow(samples_plasma_raw), "\n")

panel_dir_pla <- paste0(pla_dir, "protein-expression/", prefix_pla, "olink-explore_protein-expression_PLA-PPEA-D03_matrix_")
cardiometa_pla   <- read.csv(paste0(panel_dir_pla, "cardiometabolic.csv"))
inflammatory_pla <- read.csv(paste0(panel_dir_pla, "inflammation.csv"))
neurology_pla    <- read.csv(paste0(panel_dir_pla, "neurology.csv"))
oncology_pla     <- read.csv(paste0(panel_dir_pla, "oncology.csv"))

proteomic_raw_pla <- rbind(cardiometa_pla, inflammatory_pla, neurology_pla, oncology_pla)
rownames(proteomic_raw_pla) <- proteomic_raw_pla[[1]]
proteomic_raw_pla <- proteomic_raw_pla[, -1]
cat("PLA — samples in protein-expression matrix:", ncol(proteomic_raw_pla), "\n")

samples_plasma_full <- merge(samples_plasma_raw, metadata_prot[, c(1, 11, 13, 14)])
colnames(samples_plasma_full)[colnames(samples_plasma_full) == "case_control_other_at_baseline"] <- "condition"
samples_plasma <- samples_plasma_full %>% filter(condition %in% c("Case", "Control"))


participants_plasma <- samples_plasma %>% group_by(participant_id) %>% slice(1) %>% ungroup()
participants_plasma_ids <- unique(participants_plasma$participant_id)

cat("Plasma — N participants:", nrow(participants_plasma), "\n")
cat("Mean age:", round(mean(participants_plasma$age_at_baseline, na.rm = TRUE), 1),
    "± SD", round(sd(participants_plasma$age_at_baseline, na.rm = TRUE), 1), "\n")
cat("Sex distribution:\n"); print(table(participants_plasma$sex))
cat("Condition distribution:\n"); print(table(participants_plasma$condition))

case_control_plasma_by_time <- samples_plasma %>%
  count(visit_month, condition) %>%
  tidyr::pivot_wider(names_from = condition, values_from = n, values_fill = 0) %>%
  rename(Cases = Case, Controls = Control) %>%
  mutate(Total = Cases + Controls) %>%
  arrange(visit_month)
print(case_control_plasma_by_time)

cat("\nPLA — total samples:", sum(case_control_plasma_by_time$Total), "\n")
cat("PLA — total cases:", sum(case_control_plasma_by_time$Cases), "\n")
cat("PLA — total controls:", sum(case_control_plasma_by_time$Controls), "\n")

# ---- Cross-modality overlap ---------------------------------------------------
cat("Transcriptomic ∩ CSF:", length(intersect(participants_tx_ids, participants_csf_ids)), "\n")
cat("Transcriptomic ∩ Plasma:", length(intersect(participants_tx_ids, participants_plasma_ids)), "\n")
cat("Transcriptomic ∩ (CSF ∪ Plasma):",
    length(intersect(participants_tx_ids, union(participants_csf_ids, participants_plasma_ids))), "\n")
cat("Total unique across all three modalities:",
    length(union(participants_tx_ids, union(participants_csf_ids, participants_plasma_ids))), "\n")
