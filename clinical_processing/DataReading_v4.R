#Install required packages
library("plyr")
library("dplyr")
library("table1")
library("psych")
library("tidyr")

##1. DATA READING 

setwd(here::here("data/raw/AMP_PD_Data_v3_2022"))

participants <- read.csv("clinical_metadata_v4_2023/releases_2023_v4release_1027_amp_pd_participants.csv")
duplicates <- read.csv("clinical_metadata_v4_2023/releases_2023_v4release_1027_amp_pd_participant_wgs_duplicates.csv")
case_control <- read.csv("clinical_metadata_v4_2023/releases_2023_v4release_1027_amp_pd_case_control.csv")
global_sample <- read.csv("clinical_metadata_v4_2023/releases_2023_v4release_1027_amp_pd_global_sample_inventory.csv")

##1.1 CLINICAL

abeta <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_Biospecimen_analyses_CSF_abeta_tau_ptau.csv")
  abeta_dict <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_Biospecimen_analyses_CSF_abeta_tau_ptau_dictionary.csv")

gluco <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_Biospecimen_analyses_CSF_beta_glucocerebrosidase.csv")
  gluco_dict <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_Biospecimen_analyses_CSF_beta_glucocerebrosidase_dictionary.csv")

bio_other <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_Biospecimen_analyses_other.csv")
  bio_other_dict <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_Biospecimen_analyses_other_dictionary.csv")

soma_plasma <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_Biospecimen_analyses_SomaLogic_plasma.csv")
  soma_plasma_dict <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_Biospecimen_analyses_SomaLogic_plasma_dictionary.csv")

DaTSCAN <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_DaTSCAN_SBR.csv")
  DaTSCAN_dict <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_DaTSCAN_SBR_dictionary.csv")
DaTSCAN_visual <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_DaTSCAN_visual_interpretation.csv")
  DaTSCAN__visual_dict <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_DaTSCAN_visual_interpretation_dictionary.csv")

DTI <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_DTI.csv")
  DTI_dict <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_DTI_dictionary.csv")

enrollment <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_Enrollment.csv")
  enrollment_dict <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_Enrollment_dictionary.csv")

fam_history <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_Family_History_PD.csv")
  fam_history_dict <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_Family_History_PD_dictionary.csv")

cohort <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_LBD_Cohort_Clinical_Data.csv")
  cohort_dict <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_LBD_Cohort_Clinical_Data_dictionary.csv")

cohort_path <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_LBD_Cohort_Path_Data.csv")
  cohort_path_dict <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_LBD_Cohort_Path_Data_dictionary.csv")

caffeine <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_Caffeine_history.csv")
  caffeine_dict <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_Caffeine_history_dictionary.csv")

sleepiness <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_Epworth_Sleepiness_Scale.csv")
  sleepiness_dict <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_Epworth_Sleepiness_Scale_dictionary.csv")

sleepiness_mayo <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_REM_Sleep_Behavior_Disorder_Questionnaire_Mayo.csv")
  sleepiness_mayo_dict <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_REM_Sleep_Behavior_Disorder_Questionnaire_Mayo_dictionary.csv")
  
sleepiness_kolster <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_REM_Sleep_Behavior_Disorder_Questionnaire_Stiasny_Kolster.csv")
  sleepiness_kolster_dict <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_REM_Sleep_Behavior_Disorder_Questionnaire_Stiasny_Kolster_dictionary.csv")
  
smoke <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_Smoking_and_alcohol_history.csv")
  smoke_dict <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_Smoking_and_alcohol_history_dictionary.csv")

demo <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_Demographics.csv")
  demo_dict <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_Demographics_dictionary.csv")

UPDRS <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_UPDRS.csv")
  UPDRS_dict <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_UPDRS_dictionary.csv")
UPDRS_1 <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_MDS_UPDRS_Part_I.csv")
  UPDRS_dict_1 <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_MDS_UPDRS_Part_I_dictionary.csv")
UPDRS_2 <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_MDS_UPDRS_Part_II.csv")
  UPDRS_dict_2 <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_MDS_UPDRS_Part_II_dictionary.csv")
UPDRS_3 <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_MDS_UPDRS_Part_III.csv")
  UPDRS_dict_3 <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_MDS_UPDRS_Part_III_dictionary.csv")
UPDRS_4 <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_MDS_UPDRS_Part_IV.csv")
  UPDRS_dict_4 <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_MDS_UPDRS_Part_IV_dictionary.csv")

MMSE <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_MMSE.csv")
  MMSE_dict <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_MMSE_dictionary.csv")

MOCA <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_MOCA.csv")
  MOCA_dict <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_MOCA_dictionary.csv")
  
ADL <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_Modified_Schwab___England_ADL.csv")
  ADL_dict <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_Modified_Schwab___England_ADL_dictionary.csv")

MRI <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_MRI.csv")
  MRI_dict <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_MRI_dictionary.csv")
  
medicalhist <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_PD_Medical_History.csv")
  medicalhist_dict <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_PD_Medical_History_dictionary.csv")

PDQ_39 <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_PDQ_39.csv")
  PDQ_39_dict <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_PDQ_39_dictionary.csv")
  
UPSIT <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_UPSIT.csv")
  UPSIT_dict <- read.csv("clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_UPSIT_dictionary.csv")


