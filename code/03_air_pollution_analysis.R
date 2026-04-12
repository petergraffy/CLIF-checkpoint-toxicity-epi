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

clean_names_simple <- function(df) {
  names(df) <- names(df) %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
  df
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
  )

pm25_year <- safe_read_delim(get_config_value(config, "pm25_year_file", default = file.path(exposome_dir, "pm25_county_year.csv"))) %>%
  clean_names_simple()

no2_year <- safe_read_delim(get_config_value(config, "no2_year_file", default = file.path(exposome_dir, "no2_county_year.csv"))) %>%
  clean_names_simple()

svi_year <- safe_read_delim(get_config_value(config, "svi_file", default = file.path(exposome_dir, "svi_county_year.csv"))) %>%
  clean_names_simple()

standardize_county_year <- function(df, value_name) {
  county_col <- names(df)[str_detect(names(df), "county.*fips|fips|county_code")][1]
  year_col <- names(df)[str_detect(names(df), "^year$|calendar_year|observation_year")][1]
  value_col <- names(df)[str_detect(names(df), value_name)][1]
  
  if (any(is.na(c(county_col, year_col, value_col)))) {
    stop("Could not identify county/year/value columns for ", value_name)
  }
  
  df %>%
    transmute(
      county_code = pad_fips(.data[[county_col]]),
      icu_year = as.integer(.data[[year_col]]),
      value = as.numeric(.data[[value_col]])
    )
}

pm25_link <- standardize_county_year(pm25_year, "pm25") %>%
  rename(pm25_annual = value)

no2_link <- standardize_county_year(no2_year, "no2") %>%
  rename(no2_annual = value)

svi_link <- standardize_county_year(svi_year, "svi|rpl") %>%
  rename(svi_overall = value)

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
    n_svi_linked = sum(!is.na(svi_overall))
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
    .groups = "drop"
  )

rate_ci <- function(events, n) {
  ifelse(
    is.na(n) | n == 0,
    NA_real_,
    sqrt((events / n) * (1 - events / n) / n)
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

write_csv(analysis_exposome, file.path(out_dir, "analysis_exposome_linked.csv"))
write_csv(linkage_qc, file.path(out_dir, "linkage_qc.csv"))
write_csv(yearly_exposure_summary, file.path(out_dir, "yearly_exposure_summary.csv"))
write_csv(model_results, file.path(out_dir, "air_pollution_models.csv"))
write_csv(yearly_burden, file.path(out_dir, "yearly_burden_probable.csv"))
write_csv(quartile_phenotype, file.path(out_dir, "quartile_phenotype_rates.csv"))
write_csv(quartile_mortality_probable, file.path(out_dir, "quartile_mortality_probable.csv"))

ggsave(file.path(fig_dir, "figure1_annual_probable_rate.png"), p_yearly, width = 8, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "figure1_annual_probable_rate.pdf"), p_yearly, width = 8, height = 5)
ggsave(file.path(fig_dir, "figure2_pollution_quartile_probable_rate.png"), p_quartile_phenotype, width = 8.5, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "figure2_pollution_quartile_probable_rate.pdf"), p_quartile_phenotype, width = 8.5, height = 5.5)
ggsave(file.path(fig_dir, "figure3_pollution_quartile_mortality_probable.png"), p_quartile_mortality, width = 8.5, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "figure3_pollution_quartile_mortality_probable.pdf"), p_quartile_mortality, width = 8.5, height = 5.5)

if (exists("p_forest")) {
  ggsave(file.path(fig_dir, "figure4_air_pollution_forest.png"), p_forest, width = 8.5, height = 4.8, dpi = 300)
  ggsave(file.path(fig_dir, "figure4_air_pollution_forest.pdf"), p_forest, width = 8.5, height = 4.8)
}

print(linkage_qc)
print(model_results)
