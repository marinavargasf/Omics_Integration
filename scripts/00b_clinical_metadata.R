################################################################################
# 00b_clinical_metadata.R
# Builds the participant-level clinical metadata used by the proteomics scripts
# (07, 08, 11) and saves it to results/clinical/metadata.rds.
# The same two clinical scripts are sourced again by 12-16 to build
# `clinical_data` (visit level) and `metadata` (participant level).
################################################################################

source(here::here("clinical_processing/DataReading_v4.R"))
source(here::here("clinical_processing/data_processing.R"))

# 07, 08 and 11 select metadata columns by position: c(1, 11, 13, 14)
stopifnot(colnames(metadata)[1] == "participant_id",
          setequal(colnames(metadata)[c(11, 13, 14)],
                   c("case_control_other_at_baseline", "sex", "age_at_baseline")))

dir.create(here::here("results/clinical"), recursive = TRUE, showWarnings = FALSE)
saveRDS(metadata, here::here("results/clinical/metadata.rds"))
cat("Saved results/clinical/metadata.rds:", nrow(metadata), "participants\n")
