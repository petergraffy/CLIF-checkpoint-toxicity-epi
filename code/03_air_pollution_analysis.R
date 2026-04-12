# ============================================================
# 03_air_pollution_analysis.R
# County-level air pollution linkage and exploratory models
# for checkpoint irAE ICU epidemiology
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(lubridate)
  library(broom)
  library(scales)
})

source("utils/config.R")

analysis_dir <- get_config_value(config, "analysis_dir", default = "output/checkpoint_irae_icu_epi")
exposome_dir <- get_config_value(config, "exposome_dir", default = "exposome")
out_dir <- file.path(analysis_dir, "air_pollution")
fig_dir <- file.path(out_dir, "figures")
analysis_max_year <- 2024
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pad_fips <- function(x) {
  x %>%
    as.character() %>%
    str_extract("\\d+") %>%
    str_pad(width = 5, side = "left", pad = "0")
}

safe_read_delim <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path)
  }
  readr::read_csv(path, show_col_types = FALSE)
}

analysis_path <- file.path(analysis_dir, "analysis_dataset.csv")
if (!file.exists(analysis_path)) {
  stop(
    "analysis_dataset.csv not found in ", analysis_dir, ". ",
    "Save the analysis dataset from code/02_analysis.R before running this script."
  )
}

dat <- safe_read_delim(analysis_path) %>%
  mutate(
    county_code = pad_fips(county_code),
    icu_year = as.integer(icu_year),
    icu_month = as.integer(icu_month),
    phenotype_primary = as.integer(phenotype_primary),
    phenotype_any = as.integer(phenotype_any),
    phenotype_high = as.integer(phenotype_high),
    hospital_mortality = as.integer(hospital_mortality),
    prolonged_icu_los = as.integer(prolonged_icu_los),
    any_major_support = as.integer(any_major_support)
  ) %>%
  filter(!is.na(icu_year), icu_year <= analysis_max_year)

pm25_link <- safe_read_delim(
  get_config_value(config, "pm25_year_file", default = file.path(exposome_dir, "pm25_county_year.csv"))
) %>%
  transmute(
    county_code = pad_fips(GEOID),
    icu_year = as.integer(year),
    pm25_annual = as.numeric(pm25_mean),
    pm25_filled_2024 = as.logical(pm25_filled_2024)
  ) %>%
  filter(icu_year <= analysis_max_year)

no2_link <- safe_read_delim(
  get_config_value(config, "no2_year_file", default = file.path(exposome_dir, "no2_county_year.csv"))
) %>%
  transmute(
    county_code = pad_fips(GEOID),
    icu_year = as.integer(year),
    no2_annual = as.numeric(no2_mean)
  ) %>%
  filter(icu_year <= analysis_max_year)

svi_link <- safe_read_delim(
  get_config_value(config, "svi_file", default = file.path(exposome_dir, "svi_county_year.csv"))
) %>%
  transmute(
    county_code = pad_fips(GEOID),
    icu_year = as.integer(year),
    svi_year_release = as.integer(year_release),
    svi_overall = as.numeric(svi_overall),
    svi_theme1 = as.numeric(svi_theme1),
    svi_theme2 = as.numeric(svi_theme2),
    svi_theme3 = as.numeric(svi_theme3),
    svi_theme4 = as.numeric(svi_theme4)
  ) %>%
  filter(icu_year <= analysis_max_year)

analysis_exposome <- dat %>%
  left_join(pm25_link, by = c("county_code", "icu_year")) %>%
  left_join(no2_link, by = c("county_code", "icu_year")) %>%
  left_join(svi_link, by = c("county_code", "icu_year")) %>%
  mutate(
    pm25_q = ntile(pm25_annual, 4),
    no2_q = ntile(no2_annual, 4),
    svi_q = ntile(svi_overall, 4),
    pm25_q = factor(pm25_q, levels = 1:4, labels = c("Q1 lowest", "Q2", "Q3", "Q4 highest")),
    no2_q = factor(no2_q, levels = 1:4, labels = c("Q1 lowest", "Q2", "Q3", "Q4 highest")),
    svi_q = factor(svi_q, levels = 1:4, labels = c("Q1 lowest", "Q2", "Q3", "Q4 highest"))
  )

linkage_qc <- analysis_exposome %>%
  summarise(
    n_total = n(),
    n_county_missing = sum(is.na(county_code) | county_code == ""),
    n_pm25_linked = sum(!is.na(pm25_annual)),
    n_no2_linked = sum(!is.na(no2_annual)),
    n_svi_linked = sum(!is.na(svi_overall)),
    n_pm25_filled_2024 = sum(pm25_filled_2024 %in% TRUE, na.rm = TRUE)
  )

yearly_exposure_summary <- analysis_exposome %>%
  group_by(icu_year) %>%
  summarise(
    n = n(),
    phenotype_rate = mean(phenotype_primary, na.rm = TRUE),
    mortality_rate = mean(hospital_mortality, na.rm = TRUE),
    pm25_mean = mean(pm25_annual, na.rm = TRUE),
    no2_mean = mean(no2_annual, na.rm = TRUE),
    svi_mean = mean(svi_overall, na.rm = TRUE),
    svi_theme3_mean = mean(svi_theme3, na.rm = TRUE),
    .groups = "drop"
  )

rate_ci <- function(events, n) {
  ifelse(
    is.na(n) | n == 0,
    NA_real_,
    sqrt((events / n) * (1 - events / n) / n)
  )
}

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

fit_if_possible <- function(formula_obj, data, outcome, exposure) {
  model_dat <- data %>%
    select(all.vars(formula_obj)) %>%
    drop_na()
  
  if (nrow(model_dat) < 50) {
    return(tibble(
      outcome = outcome,
      exposure = exposure,
      term = NA_character_,
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      p.value = NA_real_,
      model_n = nrow(model_dat)
    ))
  }
  
  fit <- glm(formula_obj, data = model_dat, family = binomial())
  broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(
      outcome = outcome,
      exposure = exposure,
      model_n = nobs(fit)
    )
}

model_results <- bind_rows(
  fit_if_possible(
    phenotype_primary ~ pm25_annual + age_at_admission + factor(admit_year) + svi_overall,
    analysis_exposome,
    outcome = "probable_checkpoint_irae_icu",
    exposure = "pm25_annual"
  ),
  fit_if_possible(
    phenotype_primary ~ no2_annual + age_at_admission + factor(admit_year) + svi_overall,
    analysis_exposome,
    outcome = "probable_checkpoint_irae_icu",
    exposure = "no2_annual"
  ),
  fit_if_possible(
    hospital_mortality ~ pm25_annual + age_at_admission + phenotype_primary + factor(admit_year) + svi_overall,
    analysis_exposome,
    outcome = "hospital_mortality",
    exposure = "pm25_annual"
  ),
  fit_if_possible(
    any_major_support ~ no2_annual + age_at_admission + phenotype_primary + factor(admit_year) + svi_overall,
    analysis_exposome,
    outcome = "any_major_support",
    exposure = "no2_annual"
  )
) %>%
  filter(term %in% c("pm25_annual", "no2_annual"))

probable_dat <- analysis_exposome %>%
  filter(phenotype_primary == 1) %>%
  mutate(
    row_id = row_number(),
    remaining_life_exp = remaining_life_expectancy(age_at_admission),
    yll_if_death = remaining_life_exp
  )

probable_yll_decedents <- probable_dat %>%
  filter(hospital_mortality == 1)

pm25_reference <- probable_dat %>%
  filter(!is.na(pm25_annual)) %>%
  summarise(ref = median(pm25_annual[pm25_q == "Q1 lowest"], na.rm = TRUE)) %>%
  pull(ref)

if (length(pm25_reference) == 0 || is.na(pm25_reference) || !is.finite(pm25_reference)) {
  pm25_reference <- probable_dat %>%
    summarise(ref = quantile(pm25_annual, probs = 0.25, na.rm = TRUE)) %>%
    pull(ref)
}

pm25_yll_quartiles <- probable_dat %>%
  filter(!is.na(pm25_q)) %>%
  group_by(pm25_q) %>%
  summarise(
    n_probable_cases = n(),
    n_deaths = sum(hospital_mortality, na.rm = TRUE),
    total_yll_observed = sum(if_else(hospital_mortality == 1, yll_if_death, 0), na.rm = TRUE),
    mean_yll_per_case = total_yll_observed / n_probable_cases,
    mean_yll_per_death = if_else(n_deaths > 0, total_yll_observed / n_deaths, NA_real_),
    .groups = "drop"
  ) %>%
  arrange(pm25_q) %>%
  mutate(
    baseline_yll_per_case_q1 = mean_yll_per_case[pm25_q == "Q1 lowest"][1],
    excess_yll_per_case_vs_q1 = mean_yll_per_case - baseline_yll_per_case_q1,
    excess_total_yll_vs_q1 = excess_yll_per_case_vs_q1 * n_probable_cases
  )

pm25_mortality_probable_model <- NULL
pm25_attributable_yll_patient <- tibble()
pm25_attributable_yll_yearly <- tibble()
pm25_attributable_yll_summary <- tibble()

pm25_model_dat <- probable_dat %>%
  select(row_id, hospital_mortality, pm25_annual, age_at_admission, admit_year, svi_overall, yll_if_death) %>%
  drop_na()

if (nrow(pm25_model_dat) >= 50 && sum(pm25_model_dat$hospital_mortality, na.rm = TRUE) >= 10) {
  pm25_mortality_probable_model <- glm(
    hospital_mortality ~ pm25_annual + age_at_admission + factor(admit_year) + svi_overall,
    data = pm25_model_dat,
    family = binomial()
  )
  
  observed_prob <- predict(pm25_mortality_probable_model, newdata = pm25_model_dat, type = "response")
  counterfactual_dat <- pm25_model_dat %>%
    mutate(pm25_annual = pm25_reference)
  counterfactual_prob <- predict(pm25_mortality_probable_model, newdata = counterfactual_dat, type = "response")
  
  pm25_attributable_yll_patient <- probable_dat %>%
    select(row_id, hospitalization_id, icu_year, pm25_q) %>%
    inner_join(
      pm25_model_dat %>%
        mutate(
          predicted_mortality_observed = observed_prob,
          predicted_mortality_reference = counterfactual_prob,
          attributable_excess_death_prob = pmax(0, predicted_mortality_observed - predicted_mortality_reference),
          attributable_excess_yll = attributable_excess_death_prob * yll_if_death
        ),
      by = "row_id"
    )
  
  pm25_attributable_yll_yearly <- pm25_attributable_yll_patient %>%
    group_by(icu_year) %>%
    summarise(
      attributable_excess_yll = sum(attributable_excess_yll, na.rm = TRUE),
      attributable_excess_yll_per_case = mean(attributable_excess_yll, na.rm = TRUE),
      .groups = "drop"
    )
  
  pm25_attributable_yll_summary <- pm25_attributable_yll_patient %>%
    summarise(
      pm25_reference = pm25_reference,
      n_probable_cases_modeled = n(),
      total_attributable_excess_yll = sum(attributable_excess_yll, na.rm = TRUE),
      mean_attributable_excess_yll_per_case = mean(attributable_excess_yll, na.rm = TRUE),
      median_attributable_excess_yll_per_case = median(attributable_excess_yll, na.rm = TRUE)
    )
}

yearly_burden <- analysis_exposome %>%
  group_by(icu_year) %>%
  summarise(
    n = n(),
    cases = sum(phenotype_primary, na.rm = TRUE),
    rate = cases / n,
    se = rate_ci(cases, n),
    rate_low = pmax(0, rate - 1.96 * se),
    rate_high = pmin(1, rate + 1.96 * se),
    pm25_mean = mean(pm25_annual, na.rm = TRUE),
    .groups = "drop"
  )

line_scale_factor <- with(
  yearly_burden,
  if (all(is.na(rate)) || all(is.na(n)) || max(rate, na.rm = TRUE) <= 0) {
    1
  } else {
    max(n, na.rm = TRUE) / max(rate, na.rm = TRUE)
  }
)

trend_lines <- bind_rows(
  yearly_burden %>%
    transmute(
      icu_year,
      series = "Probable irAE burden",
      value = rate
    ),
  yearly_burden %>%
    transmute(
      icu_year,
      series = "Cancer ICU volume",
      value = n / line_scale_factor
    )
)

quartile_phenotype <- bind_rows(
  analysis_exposome %>%
    filter(!is.na(pm25_q)) %>%
    group_by(exposure_group = pm25_q) %>%
    summarise(
      pollutant = "PM2.5",
      n = n(),
      events = sum(phenotype_primary, na.rm = TRUE),
      rate = events / n,
      se = rate_ci(events, n),
      rate_low = pmax(0, rate - 1.96 * se),
      rate_high = pmin(1, rate + 1.96 * se),
      .groups = "drop"
    ),
  analysis_exposome %>%
    filter(!is.na(no2_q)) %>%
    group_by(exposure_group = no2_q) %>%
    summarise(
      pollutant = "NO2",
      n = n(),
      events = sum(phenotype_primary, na.rm = TRUE),
      rate = events / n,
      se = rate_ci(events, n),
      rate_low = pmax(0, rate - 1.96 * se),
      rate_high = pmin(1, rate + 1.96 * se),
      .groups = "drop"
    )
)

quartile_mortality_probable <- bind_rows(
  analysis_exposome %>%
    filter(phenotype_primary == 1, !is.na(pm25_q)) %>%
    group_by(exposure_group = pm25_q) %>%
    summarise(
      pollutant = "PM2.5",
      n = n(),
      events = sum(hospital_mortality, na.rm = TRUE),
      rate = events / n,
      se = rate_ci(events, n),
      rate_low = pmax(0, rate - 1.96 * se),
      rate_high = pmin(1, rate + 1.96 * se),
      .groups = "drop"
    ),
  analysis_exposome %>%
    filter(phenotype_primary == 1, !is.na(no2_q)) %>%
    group_by(exposure_group = no2_q) %>%
    summarise(
      pollutant = "NO2",
      n = n(),
      events = sum(hospital_mortality, na.rm = TRUE),
      rate = events / n,
      se = rate_ci(events, n),
      rate_low = pmax(0, rate - 1.96 * se),
      rate_high = pmin(1, rate + 1.96 * se),
      .groups = "drop"
    )
)

forest_dat <- model_results %>%
  mutate(
    label = case_when(
      outcome == "probable_checkpoint_irae_icu" & exposure == "pm25_annual" ~ "Probable irAE phenotype per 1-unit PM2.5",
      outcome == "probable_checkpoint_irae_icu" & exposure == "no2_annual" ~ "Probable irAE phenotype per 1-unit NO2",
      outcome == "hospital_mortality" & exposure == "pm25_annual" ~ "Hospital mortality per 1-unit PM2.5",
      outcome == "any_major_support" & exposure == "no2_annual" ~ "Major organ support per 1-unit NO2",
      TRUE ~ paste(outcome, exposure)
    )
  )

theme_grant <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )
}

p_yearly <- yearly_burden %>%
  ggplot(aes(x = icu_year, y = rate)) +
  geom_ribbon(aes(ymin = rate_low, ymax = rate_high), fill = "#A9D6E5", alpha = 0.35) +
  geom_line(color = "#005F73", linewidth = 1.1) +
  geom_point(size = 2.2, color = "#005F73") +
  scale_y_continuous(labels = label_percent(accuracy = 0.1)) +
  scale_x_continuous(breaks = yearly_burden$icu_year) +
  labs(
    title = "Annual burden of probable irAE-like ICU phenotype",
    subtitle = "Among all cancer ICU admissions",
    x = "ICU year",
    y = "Probable phenotype rate"
  ) +
  theme_grant()

p_trend_dual <- trend_lines %>%
  ggplot(aes(x = icu_year, y = value, color = series)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.3) +
  scale_color_manual(
    values = c(
      "Probable irAE burden" = "#AE2012",
      "Cancer ICU volume" = "#005F73"
    )
  ) +
  scale_y_continuous(
    name = "Probable irAE burden",
    labels = label_percent(accuracy = 0.1),
    sec.axis = sec_axis(
      ~ . * line_scale_factor,
      name = "Cancer ICU admissions"
    )
  ) +
  scale_x_continuous(breaks = yearly_burden$icu_year) +
  labs(
    title = "Probable irAE burden versus total cancer ICU volume",
    subtitle = "Annual trends limited to ICU admissions through 2024",
    x = "ICU year",
    color = NULL
  ) +
  theme_grant()

p_quartile_phenotype <- quartile_phenotype %>%
  ggplot(aes(x = exposure_group, y = rate, fill = pollutant)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  geom_errorbar(
    aes(ymin = rate_low, ymax = rate_high),
    position = position_dodge(width = 0.7),
    width = 0.16
  ) +
  scale_y_continuous(labels = label_percent(accuracy = 0.1)) +
  scale_fill_manual(values = c("PM2.5" = "#CA6702", "NO2" = "#BB3E03")) +
  labs(
    title = "Probable irAE-like phenotype rate by pollution quartile",
    subtitle = "Higher quartiles represent higher county-level annual exposure",
    x = "Exposure quartile",
    y = "Probable phenotype rate",
    fill = NULL
  ) +
  theme_grant()

p_quartile_mortality <- quartile_mortality_probable %>%
  ggplot(aes(x = exposure_group, y = rate, fill = pollutant)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  geom_errorbar(
    aes(ymin = rate_low, ymax = rate_high),
    position = position_dodge(width = 0.7),
    width = 0.16
  ) +
  scale_y_continuous(labels = label_percent(accuracy = 0.1)) +
  scale_fill_manual(values = c("PM2.5" = "#0A9396", "NO2" = "#005F73")) +
  labs(
    title = "Hospital mortality among probable irAE cases by pollution quartile",
    subtitle = "Exploratory unadjusted rates within the probable phenotype",
    x = "Exposure quartile",
    y = "Hospital mortality",
    fill = NULL
  ) +
  theme_grant()

if (nrow(forest_dat) > 0) {
  p_forest <- forest_dat %>%
    ggplot(aes(x = estimate, y = reorder(label, estimate))) +
    geom_vline(xintercept = 1, linetype = 2, color = "gray50") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.18, color = "#6C757D") +
    geom_point(size = 2.8, color = "#AE2012") +
    scale_x_log10() +
    labs(
      title = "Exploratory adjusted associations with air pollution",
      subtitle = "Odds ratios from county-year exposure models",
      x = "Adjusted odds ratio (log scale)",
      y = NULL
    ) +
    theme_grant()
}

p_pm25_yll_quartile <- pm25_yll_quartiles %>%
  ggplot(aes(x = pm25_q, y = total_yll_observed, fill = pm25_q)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = c("Q1 lowest" = "#94D2BD", "Q2" = "#E9D8A6", "Q3" = "#EE9B00", "Q4 highest" = "#BB3E03")) +
  labs(
    title = "Observed years of life lost by annual PM2.5 quartile",
    subtitle = "Among probable irAE-like ICU cases with hospital mortality",
    x = "PM2.5 quartile",
    y = "Observed total YLL"
  ) +
  theme_grant()

p_pm25_excess_yll_quartile <- pm25_yll_quartiles %>%
  filter(pm25_q != "Q1 lowest") %>%
  ggplot(aes(x = pm25_q, y = excess_total_yll_vs_q1, fill = pm25_q)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = c("Q2" = "#E9D8A6", "Q3" = "#EE9B00", "Q4 highest" = "#BB3E03")) +
  labs(
    title = "Excess observed YLL versus lowest PM2.5 quartile",
    subtitle = "Difference in YLL per probable case relative to Q1, scaled by quartile case counts",
    x = "PM2.5 quartile",
    y = "Excess YLL vs Q1"
  ) +
  theme_grant()

if (nrow(pm25_attributable_yll_yearly) > 0) {
  p_pm25_attr_yll_yearly <- pm25_attributable_yll_yearly %>%
    ggplot(aes(x = icu_year, y = attributable_excess_yll)) +
    geom_col(fill = "#AE2012", width = 0.7) +
    scale_x_continuous(breaks = pm25_attributable_yll_yearly$icu_year) +
    labs(
      title = "Estimated PM2.5-attributable excess YLL by year",
      subtitle = paste0("Counterfactual PM2.5 set to the Q1 median (", round(pm25_reference, 2), ") among probable irAE cases"),
      x = "ICU year",
      y = "Attributable excess YLL"
    ) +
    theme_grant()
}

write_csv(analysis_exposome, file.path(out_dir, "analysis_exposome_linked.csv"))
write_csv(linkage_qc, file.path(out_dir, "linkage_qc.csv"))
write_csv(yearly_exposure_summary, file.path(out_dir, "yearly_exposure_summary.csv"))
write_csv(model_results, file.path(out_dir, "air_pollution_models.csv"))
write_csv(yearly_burden, file.path(out_dir, "yearly_burden_probable.csv"))
write_csv(quartile_phenotype, file.path(out_dir, "quartile_phenotype_rates.csv"))
write_csv(quartile_mortality_probable, file.path(out_dir, "quartile_mortality_probable.csv"))
write_csv(pm25_yll_quartiles, file.path(out_dir, "pm25_yll_quartiles_probable.csv"))

if (nrow(pm25_attributable_yll_patient) > 0) {
  write_csv(pm25_attributable_yll_patient, file.path(out_dir, "pm25_attributable_yll_patient_level.csv"))
  write_csv(pm25_attributable_yll_yearly, file.path(out_dir, "pm25_attributable_yll_yearly.csv"))
  write_csv(pm25_attributable_yll_summary, file.path(out_dir, "pm25_attributable_yll_summary.csv"))
}

ggsave(file.path(fig_dir, "figure1_annual_probable_rate.png"), p_yearly, width = 8, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "figure1_annual_probable_rate.pdf"), p_yearly, width = 8, height = 5)
ggsave(file.path(fig_dir, "figure1b_burden_vs_cancer_icu_volume.png"), p_trend_dual, width = 8.5, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "figure1b_burden_vs_cancer_icu_volume.pdf"), p_trend_dual, width = 8.5, height = 5.5)
ggsave(file.path(fig_dir, "figure2_pollution_quartile_probable_rate.png"), p_quartile_phenotype, width = 8.5, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "figure2_pollution_quartile_probable_rate.pdf"), p_quartile_phenotype, width = 8.5, height = 5.5)
ggsave(file.path(fig_dir, "figure3_pollution_quartile_mortality_probable.png"), p_quartile_mortality, width = 8.5, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "figure3_pollution_quartile_mortality_probable.pdf"), p_quartile_mortality, width = 8.5, height = 5.5)
ggsave(file.path(fig_dir, "figure5_pm25_observed_yll_quartiles.png"), p_pm25_yll_quartile, width = 8.5, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "figure5_pm25_observed_yll_quartiles.pdf"), p_pm25_yll_quartile, width = 8.5, height = 5.5)
ggsave(file.path(fig_dir, "figure6_pm25_excess_yll_vs_q1.png"), p_pm25_excess_yll_quartile, width = 8.5, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "figure6_pm25_excess_yll_vs_q1.pdf"), p_pm25_excess_yll_quartile, width = 8.5, height = 5.5)

if (exists("p_forest")) {
  ggsave(file.path(fig_dir, "figure4_air_pollution_forest.png"), p_forest, width = 8.5, height = 4.8, dpi = 300)
  ggsave(file.path(fig_dir, "figure4_air_pollution_forest.pdf"), p_forest, width = 8.5, height = 4.8)
}

if (exists("p_pm25_attr_yll_yearly")) {
  ggsave(file.path(fig_dir, "figure7_pm25_attributable_excess_yll_yearly.png"), p_pm25_attr_yll_yearly, width = 8.5, height = 5.5, dpi = 300)
  ggsave(file.path(fig_dir, "figure7_pm25_attributable_excess_yll_yearly.pdf"), p_pm25_attr_yll_yearly, width = 8.5, height = 5.5)
}

print(linkage_qc)
print(model_results)
if (nrow(pm25_attributable_yll_summary) > 0) print(pm25_attributable_yll_summary)
