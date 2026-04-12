# ============================================================
# CLIF | Reverse Phenotype for Checkpoint-Inhibitor Toxicity ICU Cohort
# PI: Peter Graffy
#
# Strategy:
#   1) Find ICU hospitalizations
#   2) Restrict to adults with cancer present on admission
#   3) Build acute irAE-like phenotype primarily from biomarkers,
#      vitals, organ support, and steroid/rescue treatment
#   4) Keep hospital diagnosis–based organ flags as supportive only
#   5) Create possible / probable / high-confidence phenotype tiers
#
# Inputs (Parquet):
#   clif_adt.parquet
#   clif_hospitalization.parquet
#   clif_hospital_diagnosis.parquet
#   clif_labs.parquet
#   clif_vitals.parquet
#   clif_respiratory_support.parquet
#   clif_crrt_therapy.parquet
#   clif_medication_admin_continuous.parquet
#   clif_medication_admin_intermittent.parquet
#
# Outputs:
#   01_cancer_icu_broad.csv
#   02_cancer_icu_irae_features.csv
#   03_cancer_icu_checkpoint_suspected.csv
#   summary_counts.csv
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(arrow)
  library(readr)
  library(stringr)
})

source("utils/config.R")

# -----------------------------
# paths
# -----------------------------
clif_dir <- get_config_value(config, "clif_dir", required = TRUE)
out_dir  <- get_config_value(config, "output_dir", default = "output/checkpoint_irae_icu")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# helpers
# -----------------------------
safe_read_parquet <- function(path) {
  if (!file.exists(path)) return(NULL)
  arrow::read_parquet(path) %>% as_tibble()
}

parse_dt <- function(x) {
  if (inherits(x, "POSIXct")) return(x)
  suppressWarnings(ymd_hms(x, quiet = TRUE))
}

clean_string <- function(x) {
  x %>%
    as.character() %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9\\.]+", " ") %>%
    str_squish()
}

to_num <- function(x) {
  suppressWarnings(as.numeric(str_extract(as.character(x), "-?\\d+\\.?\\d*")))
}

safe_max <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  suppressWarnings(max(x, na.rm = TRUE))
}

safe_min <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  suppressWarnings(min(x, na.rm = TRUE))
}

# -----------------------------
# load data
# -----------------------------
adt   <- safe_read_parquet(file.path(clif_dir, "clif_adt.parquet"))
hosp  <- safe_read_parquet(file.path(clif_dir, "clif_hospitalization.parquet"))
hdx   <- safe_read_parquet(file.path(clif_dir, "clif_hospital_diagnosis.parquet"))
labs  <- safe_read_parquet(file.path(clif_dir, "clif_labs.parquet"))
vits  <- safe_read_parquet(file.path(clif_dir, "clif_vitals.parquet"))
rsup  <- safe_read_parquet(file.path(clif_dir, "clif_respiratory_support.parquet"))
crrt  <- safe_read_parquet(file.path(clif_dir, "clif_crrt_therapy.parquet"))
medc  <- safe_read_parquet(file.path(clif_dir, "clif_medication_admin_continuous.parquet"))
medi  <- safe_read_parquet(file.path(clif_dir, "clif_medication_admin_intermittent.parquet"))

stopifnot(!is.null(adt), !is.null(hosp), !is.null(hdx))


# -----------------------------
# standardize times / strings
# -----------------------------
adt <- adt %>%
  mutate(
    in_dttm = parse_dt(in_dttm),
    out_dttm = parse_dt(out_dttm),
    location_category = clean_string(location_category)
  )

hosp <- hosp %>%
  mutate(
    admission_dttm = parse_dt(admission_dttm),
    discharge_dttm = parse_dt(discharge_dttm)
  )

hdx <- hdx %>%
  mutate(
    diagnosis_code_clean = str_to_upper(as.character(diagnosis_code)),
    diagnosis_code_format = str_to_upper(as.character(diagnosis_code_format))
  )

if (!is.null(labs)) {
  labs <- labs %>%
    mutate(
      lab_result_dttm   = parse_dt(lab_result_dttm),
      lab_collect_dttm  = parse_dt(lab_collect_dttm),
      lab_value_num     = to_num(lab_value),
      lab_category_clean = clean_string(lab_category),
      lab_name_clean     = clean_string(lab_name)
    )
}

if (!is.null(vits)) {
  vits <- vits %>%
    mutate(
      recorded_dttm = parse_dt(recorded_dttm),
      vital_category_clean = clean_string(vital_category),
      vital_value_num = suppressWarnings(as.numeric(vital_value))
    )
}

if (!is.null(rsup)) {
  rsup <- rsup %>%
    mutate(
      recorded_dttm = parse_dt(recorded_dttm),
      device_category_clean = clean_string(device_category),
      fio2_set_num = suppressWarnings(as.numeric(fio2_set))
    )
}

if (!is.null(crrt)) {
  crrt <- crrt %>%
    mutate(
      recorded_dttm = parse_dt(recorded_dttm)
    )
}

if (!is.null(medc)) {
  medc <- medc %>%
    mutate(
      admin_dttm = parse_dt(admin_dttm),
      med_name_clean = clean_string(med_name),
      med_category_clean = clean_string(med_category)
    )
}

if (!is.null(medi)) {
  medi <- medi %>%
    mutate(
      admin_dttm = parse_dt(admin_dttm),
      med_name_clean = clean_string(med_name),
      med_category_clean = clean_string(med_category),
      mar_action_group_clean = clean_string(mar_action_group)
    )
}

# ============================================================
# 1) ICU anchor
# ============================================================
icu_hosp <- adt %>%
  filter(location_category == "icu") %>%
  group_by(hospitalization_id) %>%
  summarise(
    first_icu_in = suppressWarnings(min(in_dttm, na.rm = TRUE)),
    last_icu_out = suppressWarnings(max(out_dttm, na.rm = TRUE)),
    n_icu_segments = n(),
    .groups = "drop"
  ) %>%
  mutate(
    first_icu_in = ifelse(is.infinite(first_icu_in), NA, first_icu_in),
    last_icu_out = ifelse(is.infinite(last_icu_out), NA, last_icu_out)
  ) %>%
  mutate(
    first_icu_in = as.POSIXct(first_icu_in, origin = "1970-01-01", tz = "UTC"),
    last_icu_out = as.POSIXct(last_icu_out, origin = "1970-01-01", tz = "UTC")
  )

# ============================================================
# 2) POA cancer anchor
# ============================================================
cancer_patterns_icd10 <- c(
  "^C0", "^C1", "^C2", "^C3", "^C4", "^C5", "^C6", "^C7", "^C8", "^C9",
  "^C[1-7][0-9]",
  "^C80", "^C81", "^C82", "^C83", "^C84", "^C85", "^C86", "^C88",
  "^C90", "^C91", "^C92", "^C93", "^C94", "^C95", "^C96", "^C7A", "^C7B",
  "^V10", "^Z85"
)

cancer_patterns_icd9 <- c(
  "^14[0-9]", "^15[0-9]", "^16[0-9]", "^17[0-9]", "^18[0-9]", "^19[0-9]",
  "^20[0-9]", "^21[0-9]", "^22[0-9]", "^23[0-9]",
  "^V10"
)

poa_cancer <- hdx %>%
  filter(poa_present == 1) %>%
  group_by(hospitalization_id) %>%
  summarise(
    poa_cancer = as.integer(any(
      (diagnosis_code_format == "ICD10CM" & str_detect(diagnosis_code_clean, paste(cancer_patterns_icd10, collapse = "|"))) |
        (diagnosis_code_format == "ICD9CM"  & str_detect(diagnosis_code_clean, paste(cancer_patterns_icd9,  collapse = "|")))
    )),
    cancer_codes_poa = paste(
      sort(unique(diagnosis_code_clean[
        (diagnosis_code_format == "ICD10CM" & str_detect(diagnosis_code_clean, paste(cancer_patterns_icd10, collapse = "|"))) |
          (diagnosis_code_format == "ICD9CM"  & str_detect(diagnosis_code_clean, paste(cancer_patterns_icd9,  collapse = "|")))
      ])),
      collapse = " | "
    ),
    .groups = "drop"
  ) %>%
  filter(poa_cancer == 1)

# ============================================================
# 3) Supportive hospital diagnosis flags only
# These are not the core phenotype engine.
# ============================================================
supportive_dx_flags <- hdx %>%
  group_by(hospitalization_id) %>%
  summarise(
    dx_pneumonitis = as.integer(any(str_detect(diagnosis_code_clean, "^J70|^J84|^J96|^J98"), na.rm = TRUE)),
    dx_myocarditis = as.integer(any(str_detect(diagnosis_code_clean, "^I40|^I41|^I51\\.4|^I30|^I31"), na.rm = TRUE)),
    dx_hepatitis   = as.integer(any(str_detect(diagnosis_code_clean, "^K71|^K72|^K75|^R74"), na.rm = TRUE)),
    dx_colitis     = as.integer(any(str_detect(diagnosis_code_clean, "^K52|^K51|^A09|^R19\\.7"), na.rm = TRUE)),
    dx_nephritis   = as.integer(any(str_detect(diagnosis_code_clean, "^N10|^N12|^N14|^N17|^N05"), na.rm = TRUE)),
    dx_endocrine   = as.integer(any(str_detect(diagnosis_code_clean, "^E03|^E05|^E10|^E11|^E16|^E23|^E27"), na.rm = TRUE)),
    dx_neuro       = as.integer(any(str_detect(diagnosis_code_clean, "^G04|^G70|^G72|^G61|^G62|^R56|^G93"), na.rm = TRUE)),
    dx_hlh_like    = as.integer(any(str_detect(diagnosis_code_clean, "^D76"), na.rm = TRUE)),
    .groups = "drop"
  )

# ============================================================
# 4) Broad cancer ICU cohort
# ============================================================
cancer_icu_broad <- hosp %>%
  inner_join(icu_hosp, by = "hospitalization_id") %>%
  inner_join(poa_cancer, by = "hospitalization_id") %>%
  filter(age_at_admission >= 18)

#write_csv(cancer_icu_broad, file.path(out_dir, "01_cancer_icu_broad.csv"))

# ============================================================
# 5) Respiratory support features
# ============================================================
if (!is.null(rsup)) {
  resp_flags <- rsup %>%
    semi_join(cancer_icu_broad %>% select(hospitalization_id, first_icu_in), by = "hospitalization_id") %>%
    left_join(cancer_icu_broad %>% select(hospitalization_id, first_icu_in), by = "hospitalization_id") %>%
    filter(recorded_dttm >= first_icu_in - hours(24),
           recorded_dttm <= first_icu_in + hours(72)) %>%
    group_by(hospitalization_id) %>%
    summarise(
      any_imv   = as.integer(any(device_category_clean == "imv", na.rm = TRUE)),
      any_hfnc  = as.integer(any(device_category_clean %in% c("high flow nc", "hfnc"), na.rm = TRUE)),
      any_nippv = as.integer(any(device_category_clean %in% c("nippv", "cpap", "bipap"), na.rm = TRUE)),
      max_fio2  = safe_max(fio2_set_num),
      .groups = "drop"
    )
} else {
  resp_flags <- cancer_icu_broad %>%
    transmute(
      hospitalization_id,
      any_imv = 0L,
      any_hfnc = 0L,
      any_nippv = 0L,
      max_fio2 = NA_real_
    )
}

# ============================================================
# 6) Vasopressor features
# ============================================================
if (!is.null(medc)) {
  vaso_flags <- medc %>%
    semi_join(cancer_icu_broad %>% select(hospitalization_id, first_icu_in), by = "hospitalization_id") %>%
    left_join(cancer_icu_broad %>% select(hospitalization_id, first_icu_in), by = "hospitalization_id") %>%
    filter(admin_dttm >= first_icu_in - hours(24),
           admin_dttm <= first_icu_in + hours(72)) %>%
    mutate(
      vaso_flag =
        str_detect(med_name_clean, "norepinephrine|epinephrine|phenylephrine|vasopressin|dopamine|dobutamine|milrinone|angiotensin") |
        med_category_clean %in% c(
          "norepinephrine", "epinephrine", "phenylephrine", "vasopressin",
          "dopamine", "dobutamine", "milrinone", "angiotensin"
        )
    ) %>%
    group_by(hospitalization_id) %>%
    summarise(
      any_vasopressor = as.integer(any(vaso_flag, na.rm = TRUE)),
      .groups = "drop"
    )
} else {
  vaso_flags <- cancer_icu_broad %>%
    transmute(hospitalization_id, any_vasopressor = 0L)
}

# ============================================================
# 7) CRRT features
# ============================================================
if (!is.null(crrt)) {
  crrt_flags <- crrt %>%
    semi_join(cancer_icu_broad %>% select(hospitalization_id, first_icu_in), by = "hospitalization_id") %>%
    left_join(cancer_icu_broad %>% select(hospitalization_id, first_icu_in), by = "hospitalization_id") %>%
    filter(recorded_dttm >= first_icu_in - hours(24),
           recorded_dttm <= first_icu_in + hours(72)) %>%
    distinct(hospitalization_id) %>%
    mutate(any_crrt = 1L)
} else {
  crrt_flags <- cancer_icu_broad %>%
    transmute(hospitalization_id, any_crrt = 0L)
}

# ============================================================
# 8) Steroid / rescue treatment features
# ============================================================
if (!is.null(medi)) {
  steroid_flags <- medi %>%
    semi_join(cancer_icu_broad %>% select(hospitalization_id, first_icu_in), by = "hospitalization_id") %>%
    left_join(cancer_icu_broad %>% select(hospitalization_id, first_icu_in), by = "hospitalization_id") %>%
    filter(admin_dttm >= first_icu_in - hours(24),
           admin_dttm <= first_icu_in + hours(72)) %>%
    filter(is.na(mar_action_group_clean) | mar_action_group_clean == "administered") %>%
    mutate(
      steroid_flag =
        str_detect(med_name_clean, "methylprednisolone|prednisone|prednisolone|dexamethasone|hydrocortisone") |
        med_category_clean %in% c("methylprednisolone", "prednisone", "prednisolone", "dexamethasone", "hydrocortisone"),
      rescue_imm_flag =
        str_detect(med_name_clean, "mycophenolate|infliximab|tocilizumab|abatacept|rituximab|cyclophosphamide|immune globulin|ivig")
    ) %>%
    group_by(hospitalization_id) %>%
    summarise(
      any_steroid = as.integer(any(steroid_flag, na.rm = TRUE)),
      any_rescue_immunosuppression = as.integer(any(rescue_imm_flag, na.rm = TRUE)),
      steroid_names = paste(sort(unique(med_name_clean[steroid_flag %in% TRUE])), collapse = " | "),
      .groups = "drop"
    ) %>%
    mutate(
      steroid_names = na_if(steroid_names, "")
    )
} else {
  steroid_flags <- cancer_icu_broad %>%
    transmute(
      hospitalization_id,
      any_steroid = 0L,
      any_rescue_immunosuppression = 0L,
      steroid_names = NA_character_
    )
}

# ============================================================
# 9) Lab-based inflammatory / organ injury features
# BNP dropped from the core phenotype due to absent/poor coverage.
# ============================================================
if (!is.null(labs)) {
  labs_win <- labs %>%
    semi_join(cancer_icu_broad %>% select(hospitalization_id, first_icu_in), by = "hospitalization_id") %>%
    left_join(cancer_icu_broad %>% select(hospitalization_id, first_icu_in), by = "hospitalization_id") %>%
    filter(lab_result_dttm >= first_icu_in - hours(24),
           lab_result_dttm <= first_icu_in + hours(72))
  
  labs_flags <- labs_win %>%
    mutate(
      ferritin_flag = str_detect(lab_category_clean, "ferritin") |
        str_detect(lab_name_clean, "ferritin"),
      
      crp_flag = str_detect(lab_category_clean, "crp|c reactive|c reactive protein|hs crp") |
        str_detect(lab_name_clean, "crp|c reactive|c reactive protein|hs crp"),
      
      esr_flag = str_detect(lab_category_clean, "esr|sed rate|sedimentation") |
        str_detect(lab_name_clean, "esr|sed rate|sedimentation"),
      
      trop_flag = str_detect(lab_category_clean, "troponin") |
        str_detect(lab_name_clean, "troponin"),
      
      ast_flag = str_detect(lab_category_clean, "^ast$|sgot") |
        str_detect(lab_name_clean, "^ast$|sgot"),
      
      alt_flag = str_detect(lab_category_clean, "^alt$|sgpt") |
        str_detect(lab_name_clean, "^alt$|sgpt"),
      
      bili_flag = str_detect(lab_category_clean, "bilirubin") |
        str_detect(lab_name_clean, "bilirubin"),
      
      creat_flag = str_detect(lab_category_clean, "creatinine") |
        str_detect(lab_name_clean, "creatinine"),
      
      lactate_flag = str_detect(lab_category_clean, "lactate") |
        str_detect(lab_name_clean, "lactate"),
      
      wbc_flag = str_detect(lab_category_clean, "wbc|white blood") |
        str_detect(lab_name_clean, "wbc|white blood"),
      
      eos_flag = str_detect(lab_category_clean, "eosin") |
        str_detect(lab_name_clean, "eosin"),
      
      plt_flag = str_detect(lab_category_clean, "platelet") |
        str_detect(lab_name_clean, "platelet")
    ) %>%
    group_by(hospitalization_id) %>%
    summarise(
      ferritin_max = safe_max(lab_value_num[ferritin_flag]),
      crp_max      = safe_max(lab_value_num[crp_flag]),
      esr_max      = safe_max(lab_value_num[esr_flag]),
      troponin_max = safe_max(lab_value_num[trop_flag]),
      ast_max      = safe_max(lab_value_num[ast_flag]),
      alt_max      = safe_max(lab_value_num[alt_flag]),
      bili_max     = safe_max(lab_value_num[bili_flag]),
      creat_max    = safe_max(lab_value_num[creat_flag]),
      lactate_max  = safe_max(lab_value_num[lactate_flag]),
      wbc_max      = safe_max(lab_value_num[wbc_flag]),
      eos_max      = safe_max(lab_value_num[eos_flag]),
      plt_min      = safe_min(lab_value_num[plt_flag]),
      .groups = "drop"
    ) %>%
    mutate(
      inflam_marker_positive = as.integer(
        (!is.na(ferritin_max) & ferritin_max >= 1000) |
          (!is.na(crp_max) & crp_max >= 10) |
          (!is.na(esr_max) & esr_max >= 50)
      ),
      hepatitis_lab_positive = as.integer(
        (!is.na(ast_max) & ast_max >= 200) |
          (!is.na(alt_max) & alt_max >= 200) |
          (!is.na(bili_max) & bili_max >= 2)
      ),
      myocarditis_lab_positive = as.integer(
        (!is.na(troponin_max) & troponin_max > 0)
      ),
      renal_lab_positive = as.integer(
        (!is.na(creat_max) & creat_max >= 2)
      ),
      cytokine_storm_like = as.integer(
        (!is.na(ferritin_max) & ferritin_max >= 3000) |
          (
            (!is.na(ferritin_max) & ferritin_max >= 1000) &
              (!is.na(crp_max) & crp_max >= 10) &
              (!is.na(plt_min) & plt_min < 100)
          )
      )
    )
} else {
  labs_flags <- cancer_icu_broad %>%
    transmute(
      hospitalization_id,
      ferritin_max = NA_real_,
      crp_max = NA_real_,
      esr_max = NA_real_,
      troponin_max = NA_real_,
      ast_max = NA_real_,
      alt_max = NA_real_,
      bili_max = NA_real_,
      creat_max = NA_real_,
      lactate_max = NA_real_,
      wbc_max = NA_real_,
      eos_max = NA_real_,
      plt_min = NA_real_,
      inflam_marker_positive = 0L,
      hepatitis_lab_positive = 0L,
      myocarditis_lab_positive = 0L,
      renal_lab_positive = 0L,
      cytokine_storm_like = 0L
    )
}

# ============================================================
# 10) Vitals around ICU entry
# ============================================================
if (!is.null(vits)) {
  vital_flags <- vits %>%
    semi_join(cancer_icu_broad %>% select(hospitalization_id, first_icu_in), by = "hospitalization_id") %>%
    left_join(cancer_icu_broad %>% select(hospitalization_id, first_icu_in), by = "hospitalization_id") %>%
    filter(recorded_dttm >= first_icu_in - hours(24),
           recorded_dttm <= first_icu_in + hours(24)) %>%
    group_by(hospitalization_id) %>%
    summarise(
      temp_max = safe_max(vital_value_num[vital_category == "temp_c"]),
      hr_max   = safe_max(vital_value_num[vital_category == "heart_rate"]),
      map_min  = safe_min(vital_value_num[vital_category == "map"]),
      spo2_min = safe_min(vital_value_num[vital_category == "spo2"]),
      rr_max   = safe_max(vital_value_num[vital_category == "respiratory_rate"]),
      .groups = "drop"
    ) %>%
    mutate(
      fever_flag = as.integer(!is.na(temp_max) & temp_max >= 38.0),
      shock_vitals_flag = as.integer(!is.na(map_min) & map_min < 65),
      hypoxemia_flag = as.integer(!is.na(spo2_min) & spo2_min < 90)
    )
} else {
  vital_flags <- cancer_icu_broad %>%
    transmute(
      hospitalization_id,
      temp_max = NA_real_,
      hr_max = NA_real_,
      map_min = NA_real_,
      spo2_min = NA_real_,
      rr_max = NA_real_,
      fever_flag = 0L,
      shock_vitals_flag = 0L,
      hypoxemia_flag = 0L
    )
}

# ============================================================
# 11) Merge all features and build module-based phenotype
# ============================================================
irae_features <- cancer_icu_broad %>%
  left_join(supportive_dx_flags, by = "hospitalization_id") %>%
  left_join(resp_flags, by = "hospitalization_id") %>%
  left_join(vaso_flags, by = "hospitalization_id") %>%
  left_join(crrt_flags, by = "hospitalization_id") %>%
  left_join(steroid_flags, by = "hospitalization_id") %>%
  left_join(labs_flags, by = "hospitalization_id") %>%
  left_join(vital_flags, by = "hospitalization_id") %>%
  select(-any_of("hospitalization_joined_id")) %>%
  mutate(
    across(
      c(
        starts_with("dx_"),
        any_imv, any_hfnc, any_nippv, any_vasopressor, any_crrt,
        any_steroid, any_rescue_immunosuppression,
        inflam_marker_positive, hepatitis_lab_positive,
        myocarditis_lab_positive, renal_lab_positive,
        cytokine_storm_like, fever_flag, shock_vitals_flag,
        hypoxemia_flag
      ),
      ~ replace_na(.x, 0L)
    ),
    
    # supportive dx counts only
    organ_dx_count =
      dx_pneumonitis + dx_myocarditis + dx_hepatitis + dx_colitis +
      dx_nephritis + dx_endocrine + dx_neuro + dx_hlh_like,
    
    severe_support_count =
      any_imv + any_vasopressor + any_crrt +
      as.integer(any_hfnc == 1 & hypoxemia_flag == 1),
    
    immune_signal_count =
      inflam_marker_positive + hepatitis_lab_positive +
      myocarditis_lab_positive + renal_lab_positive +
      cytokine_storm_like + fever_flag,
    
    # biomarker/treatment-based organ modules
    pulm_irae_like = as.integer(
      (
        any_imv == 1 |
          any_hfnc == 1 |
          any_nippv == 1 |
          (!is.na(spo2_min) & spo2_min < 90) |
          (!is.na(max_fio2) & max_fio2 >= 60)
      ) &
        any_steroid == 1
    ),
    
    cardiac_irae_like = as.integer(
      (
        (!is.na(troponin_max) & troponin_max > 0)
      ) &
        (
          any_vasopressor == 1 |
            (!is.na(lactate_max) & lactate_max >= 2)
        ) &
        any_steroid == 1
    ),
    
    hepatic_irae_like = as.integer(
      (
        (!is.na(ast_max) & ast_max >= 200) |
          (!is.na(alt_max) & alt_max >= 200) |
          (!is.na(bili_max) & bili_max >= 2)
      ) &
        any_steroid == 1
    ),
    
    renal_irae_like = as.integer(
      (
        (!is.na(creat_max) & creat_max >= 2) |
          any_crrt == 1
      ) &
        any_steroid == 1
    ),
    
    hyperinflammatory_irae_like = as.integer(
      (
        (!is.na(ferritin_max) & ferritin_max >= 1000) |
          (!is.na(crp_max) & crp_max >= 10) |
          (!is.na(esr_max) & esr_max >= 50) |
          cytokine_storm_like == 1
      ) &
        (
          fever_flag == 1 |
            any_vasopressor == 1 |
            any_imv == 1
        ) &
        any_steroid == 1
    ),
    
    n_irae_modules =
      pulm_irae_like +
      cardiac_irae_like +
      hepatic_irae_like +
      renal_irae_like +
      hyperinflammatory_irae_like,
    
    possible_checkpoint_irae_icu = as.integer(
      n_irae_modules >= 1 |
        (organ_dx_count >= 1 & any_steroid == 1 & severe_support_count >= 1)
    ),
    
    probable_checkpoint_irae_icu = as.integer(
      (
        n_irae_modules >= 1 &
          (any_imv == 1 | any_vasopressor == 1 | any_crrt == 1)
      ) |
        (
          cytokine_storm_like == 1 &
            any_steroid == 1 &
            (any_imv == 1 | any_vasopressor == 1)
        )
    ),
    
    high_confidence_checkpoint_irae_icu = as.integer(
      (
        n_irae_modules >= 2 &
          (any_imv == 1 | any_vasopressor == 1 | any_crrt == 1)
      ) |
        (
          any_rescue_immunosuppression == 1 &
            (any_imv == 1 | any_vasopressor == 1 | any_crrt == 1)
        )
    )
  )

#write_csv(irae_features, file.path(out_dir, "02_cancer_icu_irae_features.csv"))

checkpoint_suspected <- irae_features %>%
  filter(possible_checkpoint_irae_icu == 1)

# write_csv(
#   checkpoint_suspected,
#   file.path(out_dir, "03_cancer_icu_checkpoint_suspected.csv")
# )

# ============================================================
# 12) Summary
# ============================================================
summary_counts <- tibble(
  n_cancer_icu_broad = nrow(cancer_icu_broad),
  n_possible_checkpoint_irae_icu = sum(irae_features$possible_checkpoint_irae_icu, na.rm = TRUE),
  n_probable_checkpoint_irae_icu = sum(irae_features$probable_checkpoint_irae_icu, na.rm = TRUE),
  n_high_confidence_checkpoint_irae_icu = sum(irae_features$high_confidence_checkpoint_irae_icu, na.rm = TRUE),
  n_pulm_irae_like = sum(irae_features$pulm_irae_like, na.rm = TRUE),
  n_cardiac_irae_like = sum(irae_features$cardiac_irae_like, na.rm = TRUE),
  n_hepatic_irae_like = sum(irae_features$hepatic_irae_like, na.rm = TRUE),
  n_renal_irae_like = sum(irae_features$renal_irae_like, na.rm = TRUE),
  n_hyperinflammatory_irae_like = sum(irae_features$hyperinflammatory_irae_like, na.rm = TRUE)
)

print(summary_counts)

#write_csv(summary_counts, file.path(out_dir, "summary_counts.csv"))










