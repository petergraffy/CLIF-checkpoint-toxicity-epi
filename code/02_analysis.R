


# ============================================================
# 00_common_functions.R
# Shared helper functions for checkpoint irAE ICU epi study
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(readr)
  library(stringr)
  library(broom)
})

`%ni%` <- Negate(`%in%`)

safe_read_csv <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  readr::read_csv(path, show_col_types = FALSE)
}

mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

summarise_binary <- function(data, var, by) {
  var_q <- rlang::enquo(var)
  by_q  <- rlang::enquo(by)
  
  if (rlang::quo_is_missing(by_q)) {
    data %>%
      summarise(
        n = dplyr::n(),
        events = sum((!!var_q) == 1, na.rm = TRUE),
        pct = 100 * events / n
      )
  } else {
    data %>%
      group_by(!!by_q) %>%
      summarise(
        n = dplyr::n(),
        events = sum((!!var_q) == 1, na.rm = TRUE),
        pct = 100 * events / n,
        .groups = "drop"
      )
  }
}

summarise_continuous <- function(data, var, by) {
  var_q <- rlang::enquo(var)
  by_q  <- rlang::enquo(by)
  
  if (rlang::quo_is_missing(by_q)) {
    data %>%
      summarise(
        n_nonmissing = sum(!is.na(!!var_q)),
        mean = mean(!!var_q, na.rm = TRUE),
        sd = sd(!!var_q, na.rm = TRUE),
        median = median(!!var_q, na.rm = TRUE),
        p25 = quantile(!!var_q, 0.25, na.rm = TRUE),
        p75 = quantile(!!var_q, 0.75, na.rm = TRUE),
        min = min(!!var_q, na.rm = TRUE),
        max = max(!!var_q, na.rm = TRUE)
      )
  } else {
    data %>%
      group_by(!!by_q) %>%
      summarise(
        n_nonmissing = sum(!is.na(!!var_q)),
        mean = mean(!!var_q, na.rm = TRUE),
        sd = sd(!!var_q, na.rm = TRUE),
        median = median(!!var_q, na.rm = TRUE),
        p25 = quantile(!!var_q, 0.25, na.rm = TRUE),
        p75 = quantile(!!var_q, 0.75, na.rm = TRUE),
        min = min(!!var_q, na.rm = TRUE),
        max = max(!!var_q, na.rm = TRUE),
        .groups = "drop"
      )
  }
}

fit_logistic_model <- function(data, formula_obj, outcome_name, exposure_name) {
  fit <- glm(formula_obj, data = data, family = binomial())
  
  broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(
      outcome = outcome_name,
      exposure = exposure_name,
      model_n = nobs(fit)
    )
}

fit_linear_model <- function(data, formula_obj, outcome_name, exposure_name) {
  fit <- lm(formula_obj, data = data)
  
  broom::tidy(fit, conf.int = TRUE) %>%
    mutate(
      outcome = outcome_name,
      exposure = exposure_name,
      model_n = nobs(fit)
    )
}

make_age_band <- function(age) {
  case_when(
    is.na(age) ~ NA_character_,
    age < 40 ~ "18-39",
    age < 50 ~ "40-49",
    age < 60 ~ "50-59",
    age < 70 ~ "60-69",
    age < 80 ~ "70-79",
    TRUE ~ "80+"
  )
}

# Approximate remaining life expectancy at age using a simple embedded table.
# This is for site-level burden estimation; pooled manuscript can replace
# with official life table if desired.
remaining_life_expectancy <- function(age) {
  case_when(
    is.na(age) ~ NA_real_,
    age < 20 ~ 60.0,
    age < 25 ~ 55.5,
    age < 30 ~ 50.8,
    age < 35 ~ 46.1,
    age < 40 ~ 41.5,
    age < 45 ~ 36.9,
    age < 50 ~ 32.4,
    age < 55 ~ 28.1,
    age < 60 ~ 23.9,
    age < 65 ~ 20.0,
    age < 70 ~ 16.4,
    age < 75 ~ 13.1,
    age < 80 ~ 10.2,
    age < 85 ~ 7.8,
    age < 90 ~ 5.9,
    TRUE ~ 4.4
  )
}

# ============================================================
# 01_make_analysis_dataset.R
# Build analysis-ready hospitalization-level dataset for
# single-site federated checkpoint irAE ICU study
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(readr)
  library(stringr)
})

# -----------------------------
# paths
# -----------------------------
phenotype_dir <- "output/checkpoint_irae_icu"
out_dir <- "output/checkpoint_irae_icu_epi"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

site_name <- "UCMC"

# -----------------------------
# read phenotype outputs
# -----------------------------
broad <- cancer_icu_broad
irae  <- irae_features

# -----------------------------
# basic cleaning
# -----------------------------
analysis_dat <- irae %>%
  mutate(
    site_name = site_name,
    admission_dttm = ymd_hms(admission_dttm, quiet = TRUE),
    discharge_dttm = ymd_hms(discharge_dttm, quiet = TRUE),
    first_icu_in   = ymd_hms(first_icu_in, quiet = TRUE),
    last_icu_out   = ymd_hms(last_icu_out, quiet = TRUE),
    admit_year = year(admission_dttm),
    icu_year   = year(first_icu_in),
    age_band = make_age_band(age_at_admission),
    
    phenotype_tier = case_when(
      high_confidence_checkpoint_irae_icu == 1 ~ "high_confidence",
      probable_checkpoint_irae_icu == 1 ~ "probable",
      possible_checkpoint_irae_icu == 1 ~ "possible",
      TRUE ~ "not_flagged"
    ),
    
    phenotype_primary = case_when(
      probable_checkpoint_irae_icu == 1 ~ 1L,
      TRUE ~ 0L
    ),
    
    phenotype_any = as.integer(possible_checkpoint_irae_icu == 1),
    phenotype_high = as.integer(high_confidence_checkpoint_irae_icu == 1),
    
    # mortality/disposition proxies from hospitalization discharge fields
    discharge_name_clean = clean_string(discharge_name),
    discharge_category_clean = clean_string(discharge_category),
    
    hospital_mortality = as.integer(
      str_detect(discharge_name_clean, "death|expired|deceased") |
        str_detect(discharge_category_clean, "death|expired|deceased")
    ),
    
    hospice_discharge = as.integer(
      str_detect(discharge_name_clean, "hospice") |
        str_detect(discharge_category_clean, "hospice")
    ),
    
    home_discharge = as.integer(
      str_detect(discharge_name_clean, "home") |
        str_detect(discharge_category_clean, "home")
    ),
    
    facility_discharge = as.integer(
      str_detect(discharge_name_clean, "rehab|skilled|snf|ltach|facility|nursing") |
        str_detect(discharge_category_clean, "rehab|skilled|snf|ltach|facility|nursing")
    ),
    
    hospital_los_days = as.numeric(difftime(discharge_dttm, admission_dttm, units = "days")),
    icu_los_days = as.numeric(difftime(last_icu_out, first_icu_in, units = "days")),
    
    prolonged_hospital_los = as.integer(!is.na(hospital_los_days) & hospital_los_days >= 14),
    prolonged_icu_los = as.integer(!is.na(icu_los_days) & icu_los_days >= 7),
    
    any_major_support = as.integer(any_imv == 1 | any_vasopressor == 1 | any_crrt == 1),
    
    cancer_type_broad = case_when(
      str_detect(cancer_codes_poa, "C3[0-9]|C4[0-1]") ~ "respiratory_thoracic",
      str_detect(cancer_codes_poa, "C1[5-9]|C2[0-1]|C26") ~ "gastrointestinal",
      str_detect(cancer_codes_poa, "C5[0-8]") ~ "breast_female_genital",
      str_detect(cancer_codes_poa, "C6[0-8]") ~ "genitourinary_male",
      str_detect(cancer_codes_poa, "C8[1-9]|C9[0-6]") ~ "hematologic",
      str_detect(cancer_codes_poa, "C4[3-9]") ~ "skin_soft_tissue",
      TRUE ~ "other_or_mixed"
    )
  ) %>%
  select(
    site_name,
    patient_id, hospitalization_id,
    admission_dttm, discharge_dttm, first_icu_in, last_icu_out,
    admit_year, icu_year,
    age_at_admission, age_band,
    admission_type_name, admission_type_category,
    discharge_name, discharge_category,
    discharge_name_clean, discharge_category_clean,
    poa_cancer, cancer_codes_poa, cancer_type_broad,
    phenotype_tier, phenotype_any, phenotype_primary, phenotype_high,
    possible_checkpoint_irae_icu, probable_checkpoint_irae_icu,
    high_confidence_checkpoint_irae_icu,
    pulm_irae_like, cardiac_irae_like, hepatic_irae_like,
    renal_irae_like, hyperinflammatory_irae_like,
    n_irae_modules,
    dx_pneumonitis, dx_myocarditis, dx_hepatitis, dx_colitis,
    dx_nephritis, dx_endocrine, dx_neuro, dx_hlh_like,
    organ_dx_count, severe_support_count, immune_signal_count,
    any_imv, any_hfnc, any_nippv, max_fio2,
    any_vasopressor, any_crrt,
    any_steroid, any_rescue_immunosuppression, steroid_names,
    ferritin_max, crp_max, esr_max, troponin_max,
    ast_max, alt_max, bili_max, creat_max, lactate_max,
    wbc_max, eos_max, plt_min,
    temp_max, hr_max, map_min, spo2_min, rr_max,
    inflam_marker_positive, hepatitis_lab_positive,
    myocarditis_lab_positive, renal_lab_positive,
    cytokine_storm_like, fever_flag, shock_vitals_flag, hypoxemia_flag,
    hospital_mortality, hospice_discharge, home_discharge,
    facility_discharge, hospital_los_days, icu_los_days,
    prolonged_hospital_los, prolonged_icu_los,
    any_major_support
  )

# -----------------------------
# denominators for burden
# -----------------------------
denominators_overall <- tibble(
  site_name = site_name,
  denominator_name = c("cancer_icu_broad"),
  n = c(nrow(analysis_dat))
)

denominators_yearly <- analysis_dat %>%
  group_by(site_name, icu_year) %>%
  summarise(
    denominator_name = "cancer_icu_broad",
    n = n(),
    .groups = "drop"
  )

# -----------------------------
# save
# -----------------------------
#write_csv(analysis_dat, file.path(out_dir, "analysis_dataset.csv"))
#write_csv(denominators_overall, file.path(out_dir, "denominators_overall.csv"))
#write_csv(denominators_yearly, file.path(out_dir, "denominators_yearly.csv"))

# -----------------------------
# quick QC
# -----------------------------
qc_summary <- tibble(
  site_name = site_name,
  n_cancer_icu_broad = nrow(analysis_dat),
  n_possible = sum(analysis_dat$possible_checkpoint_irae_icu, na.rm = TRUE),
  n_probable = sum(analysis_dat$probable_checkpoint_irae_icu, na.rm = TRUE),
  n_high_confidence = sum(analysis_dat$high_confidence_checkpoint_irae_icu, na.rm = TRUE),
  hospital_mortality_overall = mean(analysis_dat$hospital_mortality, na.rm = TRUE),
  hospital_mortality_probable = mean(analysis_dat$hospital_mortality[analysis_dat$phenotype_primary == 1], na.rm = TRUE)
)

print(qc_summary)
#write_csv(qc_summary, file.path(out_dir, "qc_summary.csv"))

# ============================================================
# 02_burden_analysis.R
# Single-site burden analysis for checkpoint irAE ICU study
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(scales)
  library(ggplot2)
})



in_dir <- "output/checkpoint_irae_icu_epi"
out_dir <- file.path(in_dir, "burden")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat <- analysis_dat
denom_yearly <- denominators_yearly

# -----------------------------
# overall burden
# -----------------------------
burden_overall <- tibble(
  site_name = unique(dat$site_name),
  denominator = "cancer_icu_broad",
  n_denominator = nrow(dat),
  n_possible = sum(dat$possible_checkpoint_irae_icu, na.rm = TRUE),
  n_probable = sum(dat$probable_checkpoint_irae_icu, na.rm = TRUE),
  n_high_confidence = sum(dat$high_confidence_checkpoint_irae_icu, na.rm = TRUE)
) %>%
  mutate(
    pct_possible = 100 * n_possible / n_denominator,
    pct_probable = 100 * n_probable / n_denominator,
    pct_high_confidence = 100 * n_high_confidence / n_denominator,
    rate_possible_per_1000 = 1000 * n_possible / n_denominator,
    rate_probable_per_1000 = 1000 * n_probable / n_denominator,
    rate_high_confidence_per_1000 = 1000 * n_high_confidence / n_denominator
  )

print(burden_overall)
#write_csv(burden_overall, file.path(out_dir, "burden_overall.csv"))

# -----------------------------
# yearly burden
# -----------------------------
cases_yearly <- dat %>%
  group_by(site_name, icu_year) %>%
  summarise(
    n_possible = sum(possible_checkpoint_irae_icu, na.rm = TRUE),
    n_probable = sum(probable_checkpoint_irae_icu, na.rm = TRUE),
    n_high_confidence = sum(high_confidence_checkpoint_irae_icu, na.rm = TRUE),
    .groups = "drop"
  )

burden_yearly <- denom_yearly %>%
  rename(n_denominator = n) %>%
  left_join(cases_yearly, by = c("site_name", "icu_year")) %>%
  mutate(
    across(c(n_possible, n_probable, n_high_confidence), ~ replace_na(.x, 0)),
    rate_possible_per_1000 = 1000 * n_possible / n_denominator,
    rate_probable_per_1000 = 1000 * n_probable / n_denominator,
    rate_high_confidence_per_1000 = 1000 * n_high_confidence / n_denominator,
    pct_possible = 100 * n_possible / n_denominator,
    pct_probable = 100 * n_probable / n_denominator,
    pct_high_confidence = 100 * n_high_confidence / n_denominator
  )

print(burden_yearly, n = 100)
#write_csv(burden_yearly, file.path(out_dir, "burden_yearly.csv"))

# -----------------------------
# phenotype composition
# -----------------------------
module_prevalence <- dat %>%
  filter(probable_checkpoint_irae_icu == 1) %>%
  summarise(
    pulm_irae_like = sum(pulm_irae_like, na.rm = TRUE),
    cardiac_irae_like = sum(cardiac_irae_like, na.rm = TRUE),
    hepatic_irae_like = sum(hepatic_irae_like, na.rm = TRUE),
    renal_irae_like = sum(renal_irae_like, na.rm = TRUE),
    hyperinflammatory_irae_like = sum(hyperinflammatory_irae_like, na.rm = TRUE)
  ) %>%
  pivot_longer(everything(), names_to = "module", values_to = "n") %>%
  mutate(
    denominator = sum(dat$probable_checkpoint_irae_icu, na.rm = TRUE),
    pct = 100 * n / denominator
  )

#write_csv(module_prevalence, file.path(out_dir, "module_prevalence_probable.csv"))

# -----------------------------
# descriptive table by phenotype tier
# -----------------------------
tier_counts <- dat %>%
  count(phenotype_tier, name = "n")

#write_csv(tier_counts, file.path(out_dir, "tier_counts.csv"))

# -----------------------------
# plots
# -----------------------------
p1 <- burden_yearly %>%
  ggplot(aes(x = icu_year, y = rate_probable_per_1000)) +
  geom_line() +
  geom_point() +
  theme_minimal(base_size = 12) +
  labs(
    title = "Annual burden of probable checkpoint irAE ICU phenotype",
    x = "ICU year",
    y = "Rate per 1,000 cancer ICU admissions"
  )

ggsave(file.path(out_dir, "burden_probable_yearly.png"), p1, width = 8, height = 5, dpi = 300)

p2 <- module_prevalence %>%
  ggplot(aes(x = reorder(module, pct), y = pct)) +
  geom_col() +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(
    title = "Module composition among probable phenotype cases",
    x = NULL,
    y = "% of probable cases"
  )

ggsave(file.path(out_dir, "module_prevalence_probable.png"), p2, width = 8, height = 5, dpi = 300)


# ============================================================
# 03_outcomes_analysis.R
# Single-site outcomes analysis for checkpoint irAE ICU study
# Primary comparison:
#   probable phenotype vs other cancer ICU admissions
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(broom)
})

in_dir <- "output/checkpoint_irae_icu_epi"
out_dir <- file.path(in_dir, "outcomes")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# define primary analysis groups
# -----------------------------
dat <- dat %>%
  mutate(
    primary_exposure = as.integer(probable_checkpoint_irae_icu == 1),
    exposure_label = if_else(primary_exposure == 1, "probable", "other_cancer_icu")
  )

# -----------------------------
# unadjusted descriptive outcomes
# -----------------------------
outcome_desc_binary <- bind_rows(
  summarise_binary(dat, hospital_mortality, exposure_label) %>% mutate(outcome = "hospital_mortality"),
  summarise_binary(dat, any_imv, exposure_label) %>% mutate(outcome = "any_imv"),
  summarise_binary(dat, any_vasopressor, exposure_label) %>% mutate(outcome = "any_vasopressor"),
  summarise_binary(dat, any_crrt, exposure_label) %>% mutate(outcome = "any_crrt"),
  summarise_binary(dat, prolonged_hospital_los, exposure_label) %>% mutate(outcome = "prolonged_hospital_los"),
  summarise_binary(dat, prolonged_icu_los, exposure_label) %>% mutate(outcome = "prolonged_icu_los"),
  summarise_binary(dat, home_discharge, exposure_label) %>% mutate(outcome = "home_discharge"),
  summarise_binary(dat, facility_discharge, exposure_label) %>% mutate(outcome = "facility_discharge"),
  summarise_binary(dat, hospice_discharge, exposure_label) %>% mutate(outcome = "hospice_discharge")
)

outcome_desc_cont <- bind_rows(
  summarise_continuous(dat, hospital_los_days, exposure_label) %>% mutate(outcome = "hospital_los_days"),
  summarise_continuous(dat, icu_los_days, exposure_label) %>% mutate(outcome = "icu_los_days"),
  summarise_continuous(dat, age_at_admission, exposure_label) %>% mutate(outcome = "age_at_admission")
)

#write_csv(outcome_desc_binary, file.path(out_dir, "outcome_descriptives_binary.csv"))
#write_csv(outcome_desc_cont, file.path(out_dir, "outcome_descriptives_continuous.csv"))

# -----------------------------
# Table 1 style characteristics
# -----------------------------
table1_binary <- bind_rows(
  summarise_binary(dat, any_imv, exposure_label) %>% mutate(variable = "any_imv"),
  summarise_binary(dat, any_vasopressor, exposure_label) %>% mutate(variable = "any_vasopressor"),
  summarise_binary(dat, any_crrt, exposure_label) %>% mutate(variable = "any_crrt"),
  summarise_binary(dat, any_steroid, exposure_label) %>% mutate(variable = "any_steroid"),
  summarise_binary(dat, any_rescue_immunosuppression, exposure_label) %>% mutate(variable = "any_rescue_immunosuppression"),
  summarise_binary(dat, pulm_irae_like, exposure_label) %>% mutate(variable = "pulm_irae_like"),
  summarise_binary(dat, cardiac_irae_like, exposure_label) %>% mutate(variable = "cardiac_irae_like"),
  summarise_binary(dat, hepatic_irae_like, exposure_label) %>% mutate(variable = "hepatic_irae_like"),
  summarise_binary(dat, renal_irae_like, exposure_label) %>% mutate(variable = "renal_irae_like"),
  summarise_binary(dat, hyperinflammatory_irae_like, exposure_label) %>% mutate(variable = "hyperinflammatory_irae_like")
)

table1_cont <- bind_rows(
  summarise_continuous(dat, age_at_admission, exposure_label) %>% mutate(variable = "age_at_admission"),
  summarise_continuous(dat, n_irae_modules, exposure_label) %>% mutate(variable = "n_irae_modules"),
  summarise_continuous(dat, ferritin_max, exposure_label) %>% mutate(variable = "ferritin_max"),
  summarise_continuous(dat, crp_max, exposure_label) %>% mutate(variable = "crp_max"),
  summarise_continuous(dat, troponin_max, exposure_label) %>% mutate(variable = "troponin_max"),
  summarise_continuous(dat, creat_max, exposure_label) %>% mutate(variable = "creat_max"),
  summarise_continuous(dat, lactate_max, exposure_label) %>% mutate(variable = "lactate_max")
)

#write_csv(table1_binary, file.path(out_dir, "table1_binary.csv"))
#write_csv(table1_cont, file.path(out_dir, "table1_continuous.csv"))

# -----------------------------
# adjusted models
# -----------------------------
# Primary adjustment set
model_dat <- dat %>%
  mutate(
    admission_type_category = fct_explicit_na(as.factor(admission_type_category), na_level = "missing"),
    cancer_type_broad = fct_explicit_na(as.factor(cancer_type_broad), na_level = "missing"),
    age_at_admission = as.numeric(age_at_admission),
    icu_year = as.numeric(icu_year)
  )

# logistic models
m_hosp_mort <- fit_logistic_model(
  model_dat,
  hospital_mortality ~ primary_exposure + age_at_admission + cancer_type_broad + admission_type_category + icu_year,
  outcome_name = "hospital_mortality",
  exposure_name = "primary_exposure"
)

m_imv <- fit_logistic_model(
  model_dat,
  any_imv ~ primary_exposure + age_at_admission + cancer_type_broad + admission_type_category + icu_year,
  outcome_name = "any_imv",
  exposure_name = "primary_exposure"
)

m_vaso <- fit_logistic_model(
  model_dat,
  any_vasopressor ~ primary_exposure + age_at_admission + cancer_type_broad + admission_type_category + icu_year,
  outcome_name = "any_vasopressor",
  exposure_name = "primary_exposure"
)

m_crrt <- fit_logistic_model(
  model_dat,
  any_crrt ~ primary_exposure + age_at_admission + cancer_type_broad + admission_type_category + icu_year,
  outcome_name = "any_crrt",
  exposure_name = "primary_exposure"
)

# linear models for LOS
m_hosp_los <- fit_linear_model(
  model_dat %>% filter(!is.na(hospital_los_days)),
  hospital_los_days ~ primary_exposure + age_at_admission + cancer_type_broad + admission_type_category + icu_year,
  outcome_name = "hospital_los_days",
  exposure_name = "primary_exposure"
)

m_icu_los <- fit_linear_model(
  model_dat %>% filter(!is.na(icu_los_days)),
  icu_los_days ~ primary_exposure + age_at_admission + cancer_type_broad + admission_type_category + icu_year,
  outcome_name = "icu_los_days",
  exposure_name = "primary_exposure"
)

all_models <- bind_rows(
  m_hosp_mort,
  m_imv,
  m_vaso,
  m_crrt,
  m_hosp_los,
  m_icu_los
)

write_csv(all_models, file.path(out_dir, "model_results_all_terms.csv"))

primary_terms <- all_models %>%
  filter(term == "primary_exposure") %>%
  select(outcome, term, estimate, conf.low, conf.high, p.value, model_n)

print(primary_terms)
write_csv(primary_terms, file.path(out_dir, "model_results_primary_exposure.csv"))

# -----------------------------
# sensitivity analyses by tier
# -----------------------------
sens_tier_desc <- dat %>%
  filter(phenotype_tier != "not_flagged") %>%
  group_by(phenotype_tier) %>%
  summarise(
    n = n(),
    hospital_mortality_pct = 100 * mean(hospital_mortality, na.rm = TRUE),
    any_imv_pct = 100 * mean(any_imv, na.rm = TRUE),
    any_vasopressor_pct = 100 * mean(any_vasopressor, na.rm = TRUE),
    any_crrt_pct = 100 * mean(any_crrt, na.rm = TRUE),
    median_hospital_los = median(hospital_los_days, na.rm = TRUE),
    median_icu_los = median(icu_los_days, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(sens_tier_desc, file.path(out_dir, "sensitivity_by_tier.csv"))

# ============================================================
# 04_yll_analysis.R
# Years of life lost analysis for checkpoint irAE ICU study
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggplot2)
})

out_dir <- file.path(in_dir, "yll")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# approximate YLL among decedents
# -----------------------------
yll_dat <- dat %>%
  filter(hospital_mortality == 1) %>%
  mutate(
    age_at_death_approx = age_at_admission,
    remaining_life_exp = remaining_life_expectancy(age_at_death_approx),
    yll = remaining_life_exp
  )

#write_csv(yll_dat, file.path(out_dir, "yll_decedents.csv"))

# -----------------------------
# summary overall and by phenotype
# -----------------------------
yll_overall <- yll_dat %>%
  summarise(
    n_deaths = n(),
    total_yll = sum(yll, na.rm = TRUE),
    mean_yll_per_death = mean(yll, na.rm = TRUE),
    median_yll_per_death = median(yll, na.rm = TRUE)
  )

yll_by_exposure <- dat %>%
  mutate(
    exposure_label = case_when(
      probable_checkpoint_irae_icu == 1 ~ "probable",
      TRUE ~ "other_cancer_icu"
    )
  ) %>%
  filter(hospital_mortality == 1) %>%
  mutate(
    remaining_life_exp = remaining_life_expectancy(age_at_admission),
    yll = remaining_life_exp
  ) %>%
  group_by(exposure_label) %>%
  summarise(
    n_deaths = n(),
    total_yll = sum(yll, na.rm = TRUE),
    mean_yll_per_death = mean(yll, na.rm = TRUE),
    median_yll_per_death = median(yll, na.rm = TRUE),
    .groups = "drop"
  )

yll_by_tier <- dat %>%
  filter(phenotype_tier != "not_flagged", hospital_mortality == 1) %>%
  mutate(
    remaining_life_exp = remaining_life_expectancy(age_at_admission),
    yll = remaining_life_exp
  ) %>%
  group_by(phenotype_tier) %>%
  summarise(
    n_deaths = n(),
    total_yll = sum(yll, na.rm = TRUE),
    mean_yll_per_death = mean(yll, na.rm = TRUE),
    median_yll_per_death = median(yll, na.rm = TRUE),
    .groups = "drop"
  )

# per-case burden
case_counts <- dat %>%
  summarise(
    n_cancer_icu = n(),
    n_possible = sum(possible_checkpoint_irae_icu, na.rm = TRUE),
    n_probable = sum(probable_checkpoint_irae_icu, na.rm = TRUE),
    n_high = sum(high_confidence_checkpoint_irae_icu, na.rm = TRUE)
  )

yll_per_case <- tibble(
  group = c("possible", "probable", "high_confidence"),
  n_cases = c(case_counts$n_possible, case_counts$n_probable, case_counts$n_high)
) %>%
  left_join(
    dat %>%
      filter(hospital_mortality == 1, phenotype_tier != "not_flagged") %>%
      mutate(
        remaining_life_exp = remaining_life_expectancy(age_at_admission),
        yll = remaining_life_exp
      ) %>%
      group_by(phenotype_tier) %>%
      summarise(total_yll = sum(yll, na.rm = TRUE), .groups = "drop") %>%
      mutate(group = phenotype_tier) %>%
      select(group, total_yll),
    by = "group"
  ) %>%
  mutate(
    total_yll = replace_na(total_yll, 0),
    yll_per_case = total_yll / n_cases
  )

write_csv(yll_overall, file.path(out_dir, "yll_overall.csv"))
write_csv(yll_by_exposure, file.path(out_dir, "yll_by_exposure.csv"))
write_csv(yll_by_tier, file.path(out_dir, "yll_by_tier.csv"))
write_csv(yll_per_case, file.path(out_dir, "yll_per_case.csv"))

print(yll_overall)
print(yll_by_exposure)
print(yll_by_tier)

# -----------------------------
# plots
# -----------------------------
if (nrow(yll_by_tier) > 0) {
  p1 <- yll_by_tier %>%
    ggplot(aes(x = phenotype_tier, y = total_yll)) +
    geom_col() +
    theme_minimal(base_size = 12) +
    labs(
      title = "Total years of life lost by phenotype tier",
      x = "Phenotype tier",
      y = "Total YLL"
    )
  
  ggsave(file.path(out_dir, "yll_total_by_tier.png"), p1, width = 7, height = 5, dpi = 300)
}

if (nrow(yll_dat) > 0) {
  p2 <- yll_dat %>%
    ggplot(aes(x = yll)) +
    geom_histogram(bins = 25) +
    theme_minimal(base_size = 12) +
    labs(
      title = "Distribution of years of life lost among decedents",
      x = "Years of life lost",
      y = "Count"
    )
  
  ggsave(file.path(out_dir, "yll_distribution.png"), p2, width = 7, height = 5, dpi = 300)
}









