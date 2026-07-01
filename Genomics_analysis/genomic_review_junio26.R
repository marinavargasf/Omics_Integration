################################################################################
# Mutation frequency by omics dataset
# Marina Vargas-Fernández
################################################################################

library(dplyr)
library(tidyr)

setwd("~/Documents/Parkinson/")

# ---- Load mutation table (one row per participant) --------------------------
# Already loaded as 'mutations' with columns:
# participant_id, has_known_GBA_mutation_in_WGS, has_known_LRRK2_mutation_in_WGS,
# has_known_SNCA_mutation_in_WGS, has_known_APOE_E4_mutation_in_WGS, has_known_PD_Mutation_in_WGS
mutations <- read.csv("clinical_metadata_v4_2023/tier2/releases_2023_v4release_1027_amp_pd_participant_mutations.csv")

# ---- Load omics participant IDs --------------------------------------------

# CSF proteomics
csf <- as.data.frame(readRDS("proteomic/filtered_proteins_metadata_CSF_SYMBOL.rds"))
csf_ids <- csf %>%
  filter(rowSums(!is.na(dplyr::select(., -participant_id, -visit_month))) > 0) %>%
  pull(participant_id) %>%
  unique()
cat("CSF participants with at least one measurement:", length(csf_ids), "\n")

# PLA proteomics
pla <- as.data.frame(readRDS("proteomic/filtered_proteins_metadata_PLA_SYMBOL.rds"))
pla_ids <- pla %>%
  filter(rowSums(!is.na(dplyr::select(., -participant_id, -visit_month))) > 0) %>%
  pull(participant_id) %>%
  unique()
cat("PLA participants with at least one measurement:", length(pla_ids), "\n")

# Transcriptomics — gene_expr must already be in environment from your RNA script
# If not, source the RNA script first or load gene_expr here
data_normalized <- as.data.frame(readRDS("rnaseq/Results_mayo/Normalization/normalized_expression.rds")) #IMPORTANTE, la que esta guardada en Results_feb25 tiene visit 0.5
rownames(data_normalized) <- sub("\\..*$", "", rownames(data_normalized))
gene_expr <- as.data.frame(t(data_normalized))
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

rna_ids <- gene_expr %>%
  filter(rowSums(!is.na(dplyr::select(., -participant_id, -visit_month, -time))) > 0) %>%
  pull(participant_id) %>%
  unique()
cat("Transcriptomics participants with at least one measurement:", length(rna_ids), "\n")

# ---- Function: mutation table for a set of IDs -----------------------------
get_mutation_counts <- function(ids, source_label, mutations_df) {
  
  df <- mutations_df %>%
    filter(participant_id %in% ids)
  
  df %>%
    dplyr::summarise(
      N_total      = n(),
      APOE_E4      = sum(has_known_APOE_E4_mutation_in_WGS == "Yes", na.rm = TRUE),
      LRRK2        = sum(has_known_LRRK2_mutation_in_WGS   == "Yes", na.rm = TRUE),
      GBA          = sum(has_known_GBA_mutation_in_WGS     == "Yes", na.rm = TRUE),
      SNCA         = sum(has_known_SNCA_mutation_in_WGS    == "Yes", na.rm = TRUE),
      Any_mutation = sum(has_known_PD_Mutation_in_WGS      == "Yes", na.rm = TRUE),
      No_mutation  = sum(has_known_APOE_E4_mutation_in_WGS == "No" &
                           has_known_LRRK2_mutation_in_WGS   == "No" &
                           has_known_GBA_mutation_in_WGS     == "No" &
                           has_known_SNCA_mutation_in_WGS    == "No", na.rm = TRUE),
      No_WGS       = sum(is.na(has_known_APOE_E4_mutation_in_WGS) |
                           is.na(has_known_LRRK2_mutation_in_WGS)   |
                           is.na(has_known_GBA_mutation_in_WGS)     |
                           is.na(has_known_SNCA_mutation_in_WGS))
    ) %>%
    pivot_longer(-N_total, names_to = "Group", values_to = "N") %>%
    mutate(
      Pct    = round(N / N_total * 100, 1),
      Source = source_label
    ) %>%
    dplyr::select(Source, Group, N, Pct)
}

# ---- Run for each omics dataset --------------------------------------------
table_csf  <- get_mutation_counts(csf_ids,  "CSF",             mutations)
table_pla  <- get_mutation_counts(pla_ids,  "PLA",             mutations)
table_rna  <- get_mutation_counts(rna_ids,  "Transcriptomics", mutations)

# ---- Combined wide table ---------------------------------------------------
combined_wide <- bind_rows(table_csf, table_pla, table_rna) %>%
  mutate(label = paste0(N, " (", Pct, "%)")) %>%
  dplyr::select(Group, Source, label) %>%
  pivot_wider(names_from = Source, values_from = label) %>%
  arrange(match(Group, c("APOE_E4", "LRRK2", "GBA", "SNCA",
                         "Any_mutation", "No_mutation", "No_WGS")))

print(combined_wide)
write.csv(combined_wide, "mutation_frequency_by_omics.csv", row.names = FALSE)
cat("Saved: mutation_frequency_by_omics.csv\n")

# ---- GBA carriers per omics dataset ----------------------------------------

get_gba_counts <- function(ids, source_label, mutations_df) {
  
  df <- mutations_df %>% filter(participant_id %in% ids)
  
  data.frame(
    Source       = source_label,
    N_total      = nrow(df),
    GBA_Yes      = sum(df$has_known_GBA_mutation_in_WGS == "Yes", na.rm = TRUE),
    GBA_No       = sum(df$has_known_GBA_mutation_in_WGS == "No",  na.rm = TRUE),
    GBA_NA       = sum(is.na(df$has_known_GBA_mutation_in_WGS))
  ) %>%
    mutate(
      Pct_Yes = round(GBA_Yes / N_total * 100, 1),
      Pct_No  = round(GBA_No  / N_total * 100, 1),
      Pct_NA  = round(GBA_NA  / N_total * 100, 1)
    )
}

gba_table <- bind_rows(
  get_gba_counts(csf_ids, "CSF",             mutations),
  get_gba_counts(pla_ids, "PLA",             mutations),
  get_gba_counts(rna_ids, "Transcriptomics", mutations)
)

print(gba_table)

