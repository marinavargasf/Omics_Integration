#Install required packages
library(qgraph)
library(bootnet) 
library(dplyr)
library(table1)
#Install required packages

#1. Data Revision

# Homogeneization of diagnosis
medicalhist$diagnosis[medicalhist$diagnosis == "Parkinson's Disease"] <- "Idiopathic PD"


# Rename of iterative visits, add 0.5 to each repetitive visit per month

add_visit_month <- function(mydata) {
  mydata$visit_month[grepl("#", mydata$visit_name)] <- 
    mydata$visit_month[grepl("#", mydata$visit_name)] + 0.5
  return(mydata)
}

abeta <- add_visit_month(abeta)
gluco <- add_visit_month(gluco)
sleepiness <- add_visit_month(sleepiness)
sleepiness_mayo <- add_visit_month(sleepiness_mayo)
sleepiness_kolster <- add_visit_month(sleepiness_kolster)
MOCA <- add_visit_month(MOCA)
UPDRS <- add_visit_month(UPDRS)
UPDRS_1 <- add_visit_month(UPDRS_1)
UPDRS_2 <- add_visit_month(UPDRS_2)
UPDRS_3 <- add_visit_month(UPDRS_3)
UPDRS_4 <- add_visit_month(UPDRS_4)
ADL <- add_visit_month(ADL)
UPSIT <- add_visit_month(UPSIT)
PDQ_39 <- add_visit_month(PDQ_39)


# Buils UPDRS summary dataset from different parts os score

F_UPDRS_1 <- UPDRS_1 %>% select(participant_id,visit_month,mds_updrs_part_i_summary_score)
F_UPDRS_2 <- UPDRS_2 %>% select(participant_id,visit_month,mds_updrs_part_ii_summary_score)
F_UPDRS_3 <- UPDRS_3 %>% select(participant_id,visit_month,mds_updrs_part_iii_summary_score, code_upd2hy_hoehn_and_yahr_stage)
F_UPDRS_4 <- UPDRS_4 %>% select(participant_id,visit_month,mds_updrs_part_iv_summary_score)

TOTAL_UPDRS <- merge(F_UPDRS_1, F_UPDRS_2, by.x = c("participant_id", "visit_month"), by.y =c("participant_id", "visit_month"), all=TRUE)

TOTAL_UPDRS <- merge(TOTAL_UPDRS, F_UPDRS_3,by.x = c("participant_id", "visit_month"), by.y =c("participant_id", "visit_month"), all=TRUE)
TOTAL_UPDRS <- merge(TOTAL_UPDRS, F_UPDRS_4,by.x = c("participant_id", "visit_month"), by.y =c("participant_id", "visit_month"), all=TRUE)

rm(F_UPDRS_1,F_UPDRS_2,F_UPDRS_3,F_UPDRS_4)

# Synthesize analytical blood test in one dataset

ptau_test <- abeta %>%filter(test_name=="p-Tau") 
colnames(ptau_test)[colnames(ptau_test)=="test_value"]<-"test_ptau"
ptau_test <- ptau_test %>% select(participant_id, visit_month, test_ptau)

tau_test <- abeta %>%filter(test_name=="Tau") 
colnames(tau_test)[colnames(tau_test)=="test_value"]<-"test_tau"
tau_test <- tau_test %>% select(participant_id, visit_month, test_tau)

abeta_test <- abeta %>%filter(test_name=="Abeta") 
colnames(abeta_test)[colnames(abeta_test)=="test_value"]<-"test_abeta"
abeta_test <- abeta_test %>% select(participant_id, visit_month, test_abeta)

table(abeta$test_name)

abeta_ptau <- merge(abeta_test, ptau_test,by.x = c("participant_id", "visit_month"), by.y =c("participant_id", "visit_month"),all = T)
abeta_multiple_test <- merge(abeta_ptau, tau_test,by.x = c("participant_id", "visit_month"), by.y =c("participant_id", "visit_month"), all=T)

abeta_multiple_test$visit_month[abeta_multiple_test$visit_month == 0.5] <- 0
rm(abeta_ptau)

# Rename cases for global family history for homogeneization

fam_history <- fam_history %>%
  filter_all(all_vars(. != "Unknown"))

fam_history <- fam_history %>%
  mutate(global_famhistory = case_when(
    biological_mother_with_pd == "Yes" | biological_father_with_pd == "Yes" | other_relative_with_pd == "Yes" ~ "Yes",
    TRUE ~ "No"
  ))

# Merge metadata 

dataframe <- merge(participants, case_control[,c(1,2,4)], all=TRUE)

dataframe <- merge(dataframe, fam_history[,c(1,8)], all=TRUE)

dataframe <- merge(dataframe, demo[,c(1,5,6,8)], all=TRUE)

metadata <- dataframe 

dim(metadata)

metadata <- metadata %>% distinct(participant_id, .keep_all = TRUE)

table(metadata$diagnosis_at_baseline)

metadata$diagnosis_at_baseline[metadata$diagnosis_at_baseline == "Parkinson's Disease"] <- "Idiopathic PD"


# Build clinical variables dataset

clinical_data <- merge(abeta_multiple_test, MOCA[,c(1,4,43)],by.x = c("participant_id", "visit_month"), by.y =c("participant_id", "visit_month"), all=TRUE)

clinical_data <- merge(clinical_data, sleepiness[,c(1,4,22)],by.x = c("participant_id", "visit_month"), by.y =c("participant_id", "visit_month"), all=TRUE)
colnames(clinical_data)[colnames(clinical_data)=="ess_summary_score"]<-"ess_sleepiness_score"

medicalhist <- medicalhist %>% filter(on_levodopa!="")
clinical_data <- merge(clinical_data, medicalhist[,c(1,4,17,18,19)],by.x = c("participant_id", "visit_month"), by.y =c("participant_id", "visit_month"), all=TRUE)

clinical_data <- merge(clinical_data, PDQ_39[,c(1,4,44:51)],by.x = c("participant_id", "visit_month"), by.y =c("participant_id", "visit_month"), all=TRUE)

clinical_data <- merge(clinical_data, TOTAL_UPDRS,by.x = c("participant_id", "visit_month"), by.y =c("participant_id", "visit_month"), all=TRUE)

clinical_data <- merge(clinical_data, ADL[,c(1,4,5)], by.x = c("participant_id", "visit_month"), by.y =c("participant_id", "visit_month"), all=TRUE)

colnames(clinical_data)[colnames(clinical_data)=="mod_schwab_england_pct_adl_score"]<-"Schwad_ADL_score"

clinical_data <- merge(clinical_data, sleepiness_kolster[,c(1,4,51)], by.x = c("participant_id", "visit_month"), by.y =c("participant_id", "visit_month"), all=TRUE)
colnames(clinical_data)[colnames(clinical_data)=="rbd_summary_score"]<-"REM_score"


clinical_data <- merge(clinical_data, UPSIT[,c(1,4,10)], by.x = c("participant_id", "visit_month"), by.y =c("participant_id", "visit_month"), all=TRUE)

clinical_data$visit_month <- as.numeric(clinical_data$visit_month)

colnames(clinical_data)

clinical_data<- clinical_data %>% select(participant_id, visit_month, test_abeta, test_ptau, test_tau, mds_updrs_part_i_summary_score,
                                          mds_updrs_part_ii_summary_score, mds_updrs_part_iii_summary_score, code_upd2hy_hoehn_and_yahr_stage
                                         ,mds_updrs_part_iv_summary_score, moca_total_score, pdq39_mobility_score, pdq39_adl_score, pdq39_emotional_score
                                         ,pdq39_stigma_score, pdq39_social_score, pdq39_cognition_score, pdq39_communication_score, pdq39_discomfort_score
                                         ,Schwad_ADL_score, REM_score, ess_sleepiness_score, upsit_total_score)


# Transform empty values to NA
clinical_data[clinical_data == ""] <- NA

metadata[metadata == ""] <- NA

# Convert all character columns to factors
metadata <- metadata %>% mutate_if(is.character, as.factor)
clinical_data <- clinical_data %>% mutate_if(is.character, as.factor)

# Fix repeated visits in same month

# Build auxiliar columns
df_preparado <- clinical_data %>%
  mutate(
    base_visit = floor(visit_month),          # Extract integer part
    is_half_visit = visit_month %% 1 == 0.5   # Verify if it is repeated visit
  )

# Identify numerical columns
numeric_cols <- names(clinical_data)[sapply(clinical_data, is.numeric)]

numeric_cols <- numeric_cols[-1]

non_numeric_cols <- names(df_preparado)[3:5]

# Grouping and calculating the mean for the numerical columns manually
df_limpio <- df_preparado %>%
  group_by(participant_id, base_visit) %>%
  dplyr::summarise(
    visit_month = first(base_visit),
    # Apply mean by rows
    across(all_of(numeric_cols), ~ mean(.x, na.rm = TRUE)),
    across(all_of(non_numeric_cols), first), # Maintain non numerical columns
    .groups = "drop"
  )

df_limpio <- as.data.frame(df_limpio)

# Rename NA for homogeneization
df_limpio <- df_limpio %>%
  mutate(across(where(is.numeric), ~ replace(., is.nan(.), NA)))

# Save results in new dataset
clinical_data <- df_limpio

# Make sure variables are coded as numeric or factors 
clinical_data$visit_month <- as.numeric(clinical_data$visit_month)
clinical_data$code_upd2hy_hoehn_and_yahr_stage <- as.factor(clinical_data$code_upd2hy_hoehn_and_yahr_stage)




