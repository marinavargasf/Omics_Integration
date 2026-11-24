# clinical_processing/

`clinical_data.R` reads and processes the AMP PD v4.0 longitudinal clinical variables. It has two parts:

1. **Data reading.** Reads the AMP PD clinical tables in `data/raw/AMP_PD_Data_v3_2022/clinical_metadata_v4_2023/`: participants, case/control status, demographics, family history, medical history and the clinical scales (MDS-UPDRS I–IV, MoCA, PDQ-39, Modified Schwab & England ADL, Epworth Sleepiness Scale, RBD questionnaire, UPSIT, …).
2. **Data processing.** Harmonises diagnoses (`"Parkinson's Disease"` → `"Idiopathic PD"`), merges the scales into one visit-level table, averages repeated visits within the same month and builds the participant-level metadata.

The file is not run on its own. It is `source()`d by `scripts/00b_clinical_metadata.R` and by the MariNET scripts 12–16. After sourcing, two objects exist in the session:

- **`clinical_data`**: one row per participant and visit (`participant_id`, `visit_month`), with the clinical scales
- **`metadata`**: one row per participant, with study, diagnosis, case/control status, demographics, family history and known PD mutations (WGS)

`00b_clinical_metadata.R` saves `metadata` to `results/clinical/metadata.rds`, which scripts 07, 08 and 11 read.

Note: the reading part calls `setwd()` on the raw clinical data folder. The pipeline scripts reset their own working directory before writing any output.
