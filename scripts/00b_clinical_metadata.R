################################################################################
# 00b_clinical_metadata.R
# Builds the participant-level clinical metadata used by the proteomics scripts
# (07, 08, 11) and saves it to results/clinical/metadata.rds.
# The same two clinical scripts are sourced again by 12-16 to build
# `clinical_data` (visit level) and `metadata` (participant level).
################################################################################

source(here::here("clinical_processing/clinical_data.R"))

# Columns used by 07, 08 and 11 (selected by name)
stopifnot(all(c("participant_id", "case_control_other_at_baseline",
                "age_at_baseline", "sex") %in% colnames(metadata)))

dir.create(here::here("results/clinical"), recursive = TRUE, showWarnings = FALSE)
saveRDS(metadata, here::here("results/clinical/metadata.rds"))
cat("Saved results/clinical/metadata.rds:", nrow(metadata), "participants\n")
