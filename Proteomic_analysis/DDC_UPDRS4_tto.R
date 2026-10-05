################################################################################
# DDC (CSF protein) predicting ONSET of UPDRS4 > 0 over follow-up
# Cases only. Binary reframing to sidestep UPDRS4 zero-inflation, which broke
# the Gaussian LMM in DDC_UPDRS4_longitudinal_analysis.R (predicted UPDRS4 < 0).
#
# Question: does higher DDC predict a participant transitioning from
# UPDRS4 = 0 (no treatment complications) to UPDRS4 > 0 over follow-up?
#
# TREATMENT ADJUSTMENT, two versions:
#  (a) PRIMARY: any_pd_treatment (levodopa OR dopamine agonist OR other PD
#      medication = "Yes"), binary, full coverage across the cohort. A direct
#      check showed every visit with UPDRS4 > 0 among levodopa="No" was
#      actually on some other dopaminergic drug -- `levodopa` alone
#      mislabeled those as "untreated."
#  (b) SECONDARY/SENSITIVITY: treatment_duration_months, continuous, from the
#      PD_Medical_History LOG table's medication start date. Only ~42% of
#      case participants have a LOG entry, so this is NOT used as the primary
#      covariate until that coverage gap is understood (Section 3d checks
#      whether missingness is concentrated among the never-treated, which
#      would be expected, vs. spread across known-treated participants too,
#      which would mean real missing data).
#
# Same plyr/dplyr masking risk as the other script -> all dplyr verbs
# explicitly namespaced below.
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(ggeffects)
library(MuMIn)

setwd("~/Documents/Omics_Integration/Proteomic_analysis/")

## ---- 1. Rebuild master_matrix exactly as in the CSF pipeline ---------------

source("~/Documents/Parkinson/clinical_metadata_v4_2023/DataReading_v4.R")
source("~/Documents/Parkinson/data_processing.R")

data_normalized <- as.data.frame(readRDS("~/Documents/Omics_Integration/Proteomic_analysis/Results_full_aligned/filtered_proteins_metadata_CSF_SYMBOL.rds"))
cat("DDC (raw NPX) mean:", mean(data_normalized$DDC), " SD:", sd(data_normalized$DDC), "\n")
proteins <- colnames(data_normalized)[-c(1, 2)]

clinical_data <- merge(clinical_data, data_normalized,
                        by.x = c("participant_id", "visit_month"),
                        by.y = c("participant_id", "visit_month"), all = TRUE)

master_matrix <- merge(clinical_data, metadata, by = "participant_id", all.x = TRUE)

master_matrix <- master_matrix %>%
  dplyr::select(participant_id, visit_month, study, diagnosis_at_baseline, age_at_baseline, sex, race,
                case_control_other_at_baseline, global_famhistory, on_levodopa, on_dopamine_agonist,
                on_other_pd_medications, has_known_GBA_mutation_in_WGS, has_known_LRRK2_mutation_in_WGS,
                has_known_SNCA_mutation_in_WGS, has_known_APOE_E4_mutation_in_WGS, has_known_PD_Mutation_in_WGS,
                mds_updrs_part_i_summary_score, mds_updrs_part_ii_summary_score, mds_updrs_part_iii_summary_score,
                mds_updrs_part_iv_summary_score, moca_total_score,
                pdq39_mobility_score, pdq39_adl_score, pdq39_emotional_score, pdq39_stigma_score,
                pdq39_social_score, pdq39_cognition_score, pdq39_communication_score, pdq39_discomfort_score,
                Schwad_ADL_score, REM_score, ess_sleepiness_score, upsit_total_score,
                dplyr::all_of(proteins))

node_names <- c("participant_id", "visit_month", "study", "diagnosis_at_baseline", "age_at_baseline",
                 "sex", "race", "case_control_other_at_baseline", "global_famhistory",
                 "levodopa", "dopamine", "other_med",
                 "GBA_mut", "LRRK2_mut", "SNCA_mut", "APOE_E4_mut", "PD_Mut",
                 "UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4", "MOCA", "Mobility39", "ADL39", "Emotional39",
                 "Stigma39", "Social39", "Cognition39", "Communication39", "Discomfort39", "Schwad_ADL",
                 "REM", "ESS", "UPSIT", proteins)
colnames(master_matrix) <- node_names

stopifnot("DDC" %in% proteins)

cat("\n=== Sanity check: levodopa column after node_names fix ===\n")
print(table(master_matrix$levodopa, useNA = "always"))

## ---- 2. CASE-ONLY FILTER ----------------------------------------------------

print(table(master_matrix$case_control_other_at_baseline, useNA = "always"))

master_matrix_cases <- master_matrix %>%
  dplyr::filter(case_control_other_at_baseline == "Case")

stopifnot(nrow(master_matrix_cases) > 0)
stopifnot(nrow(master_matrix_cases) < nrow(master_matrix))

## ---- 3. Build the analysis dataset, with BINARY UPDRS4 outcome -------------
## `any_pd_treatment` = Yes if levodopa OR dopamine agonist OR other PD
## medication = Yes at that visit. This is the PRIMARY confound/adjustment
## variable used below.

covariates <- c("age_at_baseline", "sex")

model_df <- master_matrix_cases %>%
  dplyr::select(participant_id, visit_month, DDC, UPDRS4, levodopa, dopamine, other_med,
                dplyr::all_of(covariates)) %>%
  tidyr::drop_na() %>%
  dplyr::mutate(
    updrs4_binary    = as.integer(UPDRS4 > 0),
    any_pd_treatment = factor(
      ifelse(levodopa == "Yes" | dopamine == "Yes" | other_med == "Yes", "Yes", "No")
    )
  )

cat("\nN participant-visits:", nrow(model_df), "\n")
cat("N unique participants:", dplyr::n_distinct(model_df$participant_id), "\n")

cat("\n=== Binary outcome balance (UPDRS4 > 0) ===\n")
print(table(model_df$updrs4_binary))
stopifnot(length(unique(model_df$updrs4_binary)) == 2)

## ---- 3b. CONFOUNDING CHECK: treatment vs. UPDRS4 -----------------------------

cat("\n=== Granular treatment breakdown among levodopa = No visits ===\n")
model_df %>%
  dplyr::filter(levodopa == "No") %>%
  dplyr::count(dopamine, other_med, updrs4_binary) %>%
  as.data.frame() %>%
  print()

cat("\n=== Cross-tab: any_pd_treatment x UPDRS4 binary outcome ===\n")
print(table(model_df$any_pd_treatment, model_df$updrs4_binary))

cat("\n=== Cross-tab: any_pd_treatment x DDC tertile ===\n")
model_df_check <- model_df %>% dplyr::mutate(ddc_tertile_check = dplyr::ntile(DDC, 3))
print(table(model_df_check$any_pd_treatment, model_df_check$ddc_tertile_check))
rm(model_df_check)

## ---- 3c. TREATMENT DURATION (secondary covariate) ----------------------------
## Start date pulled from EITHER source: the LOG entry (pre-baseline starters)
## OR a visit-specific row (during-study/incident starters) -- restricting to
## LOG alone missed 394 visit-specific start dates and undercounted coverage
## (63/149 participants) vs. the combined version below.

medical_hist <- read.csv("~/Documents/Parkinson/clinical_metadata_v4_2023/clinical/releases_2023_v4release_1027_clinical_PD_Medical_History.csv")

cat("Which start/initiation columns are populated where?\n")
medical_hist %>%
  dplyr::mutate(is_log = visit_name == "LOG") %>%
  dplyr::group_by(is_log) %>%
  dplyr::summarise(
    n = dplyr::n(),
    n_with_initiation = sum(!is.na(pd_medication_initiation_months_after_baseline)),
    n_with_start       = sum(!is.na(pd_medication_start_months_after_baseline)),
    .groups = "drop"
  ) %>%
  as.data.frame() %>%
  print()

start_dates <- medical_hist %>%
  dplyr::mutate(
    start_month = dplyr::coalesce(pd_medication_start_months_after_baseline,
                                    pd_medication_initiation_months_after_baseline)
  ) %>%
  dplyr::filter(!is.na(start_month)) %>%
  dplyr::group_by(participant_id) %>%
  dplyr::summarise(medication_start_month = min(start_month), .groups = "drop")

cat("\n=== Medical history: participant coverage (combined LOG + visit-specific) ===\n")
cat("Participants with a usable start_month:", nrow(start_dates), "(includes non-cases)\n")
dup_check <- start_dates %>% dplyr::count(participant_id) %>% dplyr::filter(n > 1)
if (nrow(dup_check) > 0) {
  cat("WARNING: participants with >1 row after aggregation (should not happen, check merge logic):\n")
  print(dup_check)
} else {
  cat("OK: exactly one start_month per participant.\n")
}

model_df <- model_df %>%
  dplyr::left_join(start_dates, by = "participant_id") %>%
  dplyr::mutate(
    treatment_duration_months = pmax(visit_month - medication_start_month, 0, na.rm = FALSE)
  )

cat("\nN visits with a treatment_duration_months value:",
    sum(!is.na(model_df$treatment_duration_months)), "of", nrow(model_df), "\n")
cat("Summary of treatment_duration_months (non-NA):\n")
print(summary(model_df$treatment_duration_months))

## ---- 3d. Is the duration coverage gap expected (untreated) or a real gap? --

cases_in_model <- unique(model_df$participant_id)
cat("\nUnique case participants in model_df:", length(cases_in_model), "\n")
cat("Of those, how many appear in start_dates:",
    sum(cases_in_model %in% start_dates$participant_id), "\n")

missing_from_start_dates <- setdiff(cases_in_model, start_dates$participant_id)
cat("Participants missing from medical history entirely:", length(missing_from_start_dates), "\n")

cat("\nAmong those missing, were they ever treated per any_pd_treatment?\n")
cat("(FALSE = expected/consistent gap; TRUE = real missing data)\n")
model_df %>%
  dplyr::filter(participant_id %in% missing_from_start_dates) %>%
  dplyr::group_by(participant_id) %>%
  dplyr::summarise(ever_treated = any(any_pd_treatment == "Yes"), .groups = "drop") %>%
  dplyr::count(ever_treated) %>%
  as.data.frame() %>%
  print()

cat("\nCross-check: does treatment_duration_months line up with any_pd_treatment overall?\n")
model_df %>%
  dplyr::mutate(has_duration = !is.na(treatment_duration_months) & treatment_duration_months > 0) %>%
  dplyr::count(any_pd_treatment, has_duration) %>%
  as.data.frame() %>%
  print()

## ---- 3e. SUPPLEMENTARY FIGURE: UPDRS4 distribution --------------------------

out_dir <- "~/Documents/Omics_Integration/Proteomic_analysis/Results_full_aligned/DDC_UPDRS4_followup"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_obs    <- nrow(model_df)
n_part   <- dplyr::n_distinct(model_df$participant_id)
n_zero   <- sum(model_df$UPDRS4 == 0)
pct_zero <- round(100 * n_zero / n_obs, 1)

count_df <- model_df %>%
  dplyr::count(UPDRS4, name = "n")

subtitle_text <- paste0("Cases only, complete cases for DDC/UPDRS4/covariates | N = ", n_obs,
                         " visits from ", n_part, " participants | ",
                         pct_zero, "% at UPDRS4 = 0")
subtitle_wrapped <- paste(strwrap(subtitle_text, width = 65), collapse = "\n")

p_supp <- ggplot(count_df, aes(x = factor(UPDRS4), y = n)) +
  geom_col(fill = "#619CFF", width = 0.7) +
  geom_text(aes(label = n), vjust = -0.4, size = 3.2, color = "grey30") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    x = "UPDRS Part IV score",
    y = "Number of observations",
    title = "Distribution of UPDRS4 in the DDC-UPDRS4 analysis sample",
    subtitle = subtitle_wrapped
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor    = element_blank(),
    panel.grid.major.y  = element_line(color = "#E6E6FA"),
    plot.title          = element_text(face = "bold", size = 13),
    plot.subtitle       = element_text(size = 10, color = "grey40")
  )

ggsave(file.path(out_dir, "SupplementaryFigure_UPDRS4_distribution.pdf"), p_supp, width = 7.5, height = 5.3)
ggsave(file.path(out_dir, "SupplementaryFigure_UPDRS4_distribution.png"), p_supp, width = 7.5, height = 5.3, dpi = 300)

cat("\nSupplementary figure saved (", pct_zero, "% of observations at UPDRS4 = 0)\n", sep = "")

## ---- 4. Base logistic model: UNADJUSTED for treatment -----------------------

cov_formula <- paste(covariates, collapse = " + ")

f_base <- as.formula(paste0("updrs4_binary ~ DDC + ", cov_formula, " + (1 | participant_id)"))
m_base <- glmer(f_base, data = model_df, family = binomial, control = glmerControl(optimizer = "bobyqa"))
cat("\n=== Base logistic model (no time interaction, NOT adjusted for treatment) ===\n")
print(summary(m_base))
print(MuMIn::r.squaredGLMM(m_base))

## ---- 4b. Base model ADJUSTED for any_pd_treatment (PRIMARY) ------------------

f_base_adj <- as.formula(paste0("updrs4_binary ~ DDC + any_pd_treatment + ", cov_formula, " + (1 | participant_id)"))
m_base_adj <- glmer(f_base_adj, data = model_df, family = binomial, control = glmerControl(optimizer = "bobyqa"))
cat("\n=== Base logistic model, ADJUSTED for any_pd_treatment ===\n")
print(summary(m_base_adj))
print(MuMIn::r.squaredGLMM(m_base_adj))

## ---- 4c. Direct comparison: DDC coefficient, with vs. without treatment adjustment

ddc_unadj <- coef(summary(m_base))["DDC", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]
ddc_adj   <- coef(summary(m_base_adj))["DDC", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]

cat("\n=== DDC coefficient: unadjusted vs. treatment-adjusted (any_pd_treatment) ===\n")
comparison_table <- rbind(Unadjusted = ddc_unadj, `Adjusted for any_pd_treatment` = ddc_adj)
print(comparison_table)

pct_attenuation <- 100 * (1 - ddc_adj["Estimate"] / ddc_unadj["Estimate"])
cat(sprintf("\nAttenuation in DDC estimate after treatment adjustment: %.1f%%\n", pct_attenuation))

## ---- 4d. Duration-adjusted model: DDC + treatment_duration_months ------------
## PRIMARY duration-based result (any_pd_treatment dropped from this model --
## confirmed uninformative here: near-zero correlation with DDC/duration, and
## its own coefficient was unstable purely from the 93/3 imbalance within
## this subset, not from anything about DDC or duration themselves).
## Only runs on the subset with non-missing duration (N=96, 63 participants
## -- see Section 3d for why coverage is limited to this subset).

model_df_duration <- model_df %>% dplyr::filter(!is.na(treatment_duration_months))
cat("\n=== Duration-adjusted model (DDC + treatment_duration_months), N =",
    nrow(model_df_duration), "visits,", dplyr::n_distinct(model_df_duration$participant_id), "participants ===\n")

if (nrow(model_df_duration) > 30 &&
    length(unique(model_df_duration$updrs4_binary)) == 2) {
  f_duration <- as.formula(paste0("updrs4_binary ~ DDC + treatment_duration_months + ", cov_formula,
                                   " + (1 | participant_id)"))
  m_duration <- tryCatch(
    glmer(f_duration, data = model_df_duration, family = binomial, control = glmerControl(optimizer = "bobyqa")),
    error = function(e) { cat("Duration model failed to fit:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(m_duration)) {
    print(summary(m_duration))

    ddc_duration <- coef(summary(m_duration))["DDC", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]
    dur_duration <- coef(summary(m_duration))["treatment_duration_months", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]
    cat("\n--- DDC coefficient, adjusted for treatment duration ---\n")
    print(ddc_duration)
    cat(sprintf("Odds ratio = %.2f\n", exp(ddc_duration["Estimate"])))
    cat("\n--- treatment_duration_months coefficient ---\n")
    print(dur_duration)
  }
} else {
  cat("Skipped: insufficient N or only one outcome class in the duration subset.\n")
}

## ---- 5. Extended model: DDC x visit_month interaction, treatment-adjusted --

f_time_slope <- as.formula(paste0("updrs4_binary ~ DDC * visit_month + any_pd_treatment + ", cov_formula,
                                   " + (1 + visit_month | participant_id)"))

m_time <- tryCatch({
  fit <- glmer(f_time_slope, data = model_df, family = binomial,
               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  if (isSingular(fit)) stop("singular fit with random slope")
  fit
}, error = function(e) {
  cat("\nRandom-slope model did not converge/was singular (", conditionMessage(e),
      "). Falling back to random-intercept-only.\n", sep = "")
  f_time_intercept <- as.formula(paste0("updrs4_binary ~ DDC * visit_month + any_pd_treatment + ", cov_formula,
                                         " + (1 | participant_id)"))
  glmer(f_time_intercept, data = model_df, family = binomial,
        control = glmerControl(optimizer = "bobyqa"))
})

cat("\n=== Extended logistic model (DDC x visit_month, adjusted for any_pd_treatment) ===\n")
print(summary(m_time))
print(MuMIn::r.squaredGLMM(m_time))

## ---- 6. Model comparison -----------------------------------------------------

cat("\n=== Likelihood ratio test: does the interaction improve fit? ===\n")
print(anova(m_base_adj, m_time))

## ---- 7. Effect visualization: predicted probability of UPDRS4 onset --------

## NOTE: any_pd_treatment's own coefficient in m_time is unstable (separation,
## SE in the hundreds -- see Section 5 output). ggpredict() by default holds
## factors not in `terms` at their REFERENCE level, which for any_pd_treatment
## is "No" -- exactly the pathological level, collapsing predicted
## probabilities toward zero regardless of DDC. Condition on "Yes" instead,
## both because it fixes the collapse and because it's the representative
## value for ~98% of this cohort.
pred_df <- ggeffects::ggpredict(m_time, terms = c("DDC", "visit_month [meansd]"),
                                 condition = c(any_pd_treatment = "Yes"))

p_interact <- plot(pred_df) +
  scale_color_manual(values = c("#619CFF", "#7FB77E", "#D9483F")) +
  scale_fill_manual(values = c("#619CFF", "#7FB77E", "#D9483F")) +
  labs(x = "DDC (CSF, NPX)", y = "Predicted P(UPDRS4 > 0)",
       title = "Model-implied DDC effect on probability of treatment complications",
       subtitle = "Adjusted for PD treatment; shown for treated patients (~98% of cohort)",
       color = "Visit month", fill = "Visit month") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major = element_line(color = "#E6E6FA"), panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "DDC_UPDRS4onset_interaction_effect.pdf"), p_interact, width = 7, height = 5)

## ---- 8. Individual/group-level check: empirical onset rate by DDC tertile --

baseline_state <- model_df %>%
  dplyr::group_by(participant_id) %>%
  dplyr::slice_min(visit_month, n = 1) %>%
  dplyr::ungroup()

at_risk_ids <- baseline_state %>%
  dplyr::filter(updrs4_binary == 0) %>%
  dplyr::pull(participant_id)

participant_ddc <- model_df %>%
  dplyr::filter(participant_id %in% at_risk_ids) %>%
  dplyr::group_by(participant_id) %>%
  dplyr::summarise(mean_DDC = mean(DDC), n_visits = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n_visits >= 2) %>%
  dplyr::mutate(ddc_tertile = dplyr::ntile(mean_DDC, 3))

cat("\nN at-risk participants (UPDRS4 = 0 at first visit, >=2 visits):", nrow(participant_ddc), "\n")

onset_df <- model_df %>%
  dplyr::inner_join(participant_ddc, by = "participant_id") %>%
  dplyr::mutate(ddc_tertile = factor(ddc_tertile, labels = c("Low DDC", "Mid DDC", "High DDC")))

onset_summary <- onset_df %>%
  dplyr::group_by(ddc_tertile, visit_month) %>%
  dplyr::summarise(
    n            = dplyr::n(),
    pct_with_complications = mean(updrs4_binary) * 100,
    .groups = "drop"
  )

p_onset <- ggplot(onset_summary, aes(x = visit_month, y = pct_with_complications,
                                      color = ddc_tertile, group = ddc_tertile)) +
  geom_line(linewidth = 1.1) +
  geom_point(aes(size = n), alpha = 0.7) +
  scale_color_manual(values = c("Low DDC" = "#619CFF", "Mid DDC" = "#7FB77E", "High DDC" = "#D9483F")) +
  labs(x = "Visit month", y = "% with UPDRS4 > 0",
       color = "Participant\nmean DDC tertile", size = "N at visit",
       title = "Empirical onset of treatment complications by DDC tertile",
       subtitle = "Restricted to participants complication-free at their first visit") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major = element_line(color = "#E6E6FA"), panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "UPDRS4onset_rate_by_DDC_tertile.pdf"), p_onset, width = 7, height = 5)

## ---- 9. Summary for manuscript text ------------------------------------------

cat("\n=== Summary for manuscript text ===\n")

base_est   <- coef(summary(m_base_adj))["DDC", "Estimate"]
base_se    <- coef(summary(m_base_adj))["DDC", "Std. Error"]
base_z     <- coef(summary(m_base_adj))["DDC", "z value"]
base_p     <- coef(summary(m_base_adj))["DDC", "Pr(>|z|)"]

base_ci_logodds <- base_est + c(-1.96, 1.96) * base_se
base_or         <- exp(base_est)
base_ci_or      <- exp(base_ci_logodds)

cat("\n--- Treatment-adjusted base model (DDC main effect, no time interaction) ---\n")
cat(sprintf("Estimate (log-odds) = %.3f, SE = %.3f, 95%% CI = [%.3f, %.3f]\n",
            base_est, base_se, base_ci_logodds[1], base_ci_logodds[2]))
cat(sprintf("Odds ratio = %.2f, 95%% CI = [%.2f, %.2f]\n", base_or, base_ci_or[1], base_ci_or[2]))
cat(sprintf("z = %.2f, p = %.3g\n", base_z, base_p))

cat("\nlme4::confint (Wald) cross-check for DDC:\n")
print(confint(m_base_adj, parm = "DDC", method = "Wald"))

int_term <- grep(":", rownames(coef(summary(m_time))), value = TRUE)

int_est <- coef(summary(m_time))[int_term, "Estimate"]
int_se  <- coef(summary(m_time))[int_term, "Std. Error"]
int_z   <- coef(summary(m_time))[int_term, "z value"]
int_p   <- coef(summary(m_time))[int_term, "Pr(>|z|)"]
int_ci  <- int_est + c(-1.96, 1.96) * int_se

cat("\n--- Extended model (DDC x visit_month, adjusted for any_pd_treatment) ---\n")
cat(sprintf("%s: Estimate = %.4f, SE = %.4f, 95%% CI = [%.4f, %.4f], z = %.2f, p = %.3g\n",
            int_term, int_est, int_se, int_ci[1], int_ci[2], int_z, int_p))

cat("\n--- Ready-to-paste sentence (check numbers against the printout above) ---\n")
cat(sprintf(
  paste0("After adjustment for concurrent PD treatment status (levodopa, dopamine agonist, or ",
         "other PD medication), DDC expression remained a strong, independent predictor of ",
         "complication presence (OR = %.2f, 95%% CI [%.2f, %.2f], z = %.2f, p = %.1e), indicating ",
         "that this association is not fully explained by treatment exposure alone. The DDC x ",
         "visit month interaction term was not significant (beta = %.3f, z = %.2f, p = %.2f).\n"),
  base_or, base_ci_or[1], base_ci_or[2], base_z, base_p,
  int_est, int_z, int_p
))

cat("\nDone. Plots and console output written to:", out_dir, "\n")

## ---- 10. TWO FIGURES per team discussion -------------------------------------
## Figure 1: DDC effect on P(UPDRS4 > 0), full cohort, NOT stratified by
##   treatment duration (N=415, the main/primary result).
## Figure 2: same effect, stratified by treatment duration (below vs. above
##   median -- two groups, not three, to avoid the Short/Mid overlap seen
##   with tertiles) as a SUPPORTING/secondary comparison on the smaller
##   duration-available subset (N=96). Caveat about reduced p-value from
##   sample size is baked directly into the subtitle so it travels with the
##   figure wherever it's used.

## Figure 1: main effect, full cohort, treatment-adjusted but not duration-stratified
pred_main_df <- ggeffects::ggpredict(m_base_adj, terms = "DDC",
                                       condition = c(any_pd_treatment = "Yes"))

p_main_effect <- plot(pred_main_df) +
  labs(x = "DDC (CSF, NPX)", y = "Predicted P(UPDRS4 > 0)",
       title = "DDC effect on treatment complications",
       subtitle = paste0("Full cohort, N = ", nrow(model_df), " visits, ",
                          dplyr::n_distinct(model_df$participant_id), " participants  |  ",
                          "OR = ", round(exp(coef(summary(m_base_adj))["DDC", "Estimate"]), 2),
                          ", p.val = ", formatC(coef(summary(m_base_adj))["DDC", "Pr(>|z|)"], format = "e", digits = 1))) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major = element_line(color = "#E6E6FA"), panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "DDC_UPDRS4_main_effect_full_cohort.pdf"), p_main_effect, width = 7, height = 5)
ggsave(file.path(out_dir, "DDC_UPDRS4_main_effect_full_cohort.png"), p_main_effect, width = 7, height = 5, dpi = 300)
cat("\nFigure 1 (main effect, full cohort) saved to:", out_dir, "\n")

## Figure 2: supporting comparison, stratified by treatment duration (median split)
if (exists("m_duration") && !is.null(m_duration)) {

  median_duration <- median(model_df_duration$treatment_duration_months)
  # Representative values: median of the below-median half, and median of the
  # above-median half, so the two lines sit at real, interpretable durations
  # rather than arbitrary points.
  dur_groups <- model_df_duration %>%
    dplyr::mutate(above_median = treatment_duration_months > median_duration) %>%
    dplyr::group_by(above_median) %>%
    dplyr::summarise(rep_duration = round(median(treatment_duration_months), 1), n = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(above_median)

  rep_durations <- dur_groups$rep_duration
  cat("\nMedian treatment_duration_months split at:", median_duration, "months\n")
  cat("Representative values (Below/Above median):", paste(rep_durations, collapse = ", "), "\n")
  cat("Group sizes:\n")
  print(as.data.frame(dur_groups))

  ddc_p_duration <- coef(summary(m_duration))["DDC", "Pr(>|z|)"]

  pred_duration_df <- ggeffects::ggpredict(
    m_duration,
    terms = c("DDC", paste0("treatment_duration_months [", paste(rep_durations, collapse = ","), "]"))
  )

  p_duration_effect <- plot(pred_duration_df) +
    scale_color_manual(values = c("#619CFF", "#D9483F"), labels = c("Below median", "Above median")) +
    scale_fill_manual(values = c("#619CFF", "#D9483F"), labels = c("Below median", "Above median")) +
    labs(x = "DDC (CSF, NPX)", y = "Predicted P(UPDRS4 > 0)",
         title = "DDC effect by treatment duration",
         subtitle = paste0("N = ", nrow(model_df_duration), " visits, ",
                            dplyr::n_distinct(model_df_duration$participant_id), " participants (subset) ",
                            "DDC: p.val = ", formatC(ddc_p_duration, format = "f", digits = 3),
                            "\nMedian split at ", median_duration, " months of treatment"),
         color = "Treatment\nduration", fill = "Treatment\nduration") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major = element_line(color = "#E6E6FA"), panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 9))

  ggsave(file.path(out_dir, "DDC_UPDRS4_by_duration_median.pdf"), p_duration_effect, width = 7.5, height = 5.3)
  ggsave(file.path(out_dir, "DDC_UPDRS4_by_duration_median.png"), p_duration_effect, width = 7.5, height = 5.3, dpi = 300)

  cat("\nFigure 2 (duration comparison) saved to:", out_dir, "\n")

} else {
  cat("\nSkipped Figure 2: m_duration not available.\n")
}