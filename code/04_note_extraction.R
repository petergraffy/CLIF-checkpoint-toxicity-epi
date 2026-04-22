# ============================================================
# 04_note_extraction.R
# Hybrid note extraction workflow for ICI irAE chart review
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(lubridate)
})

find_first_existing <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) {
    return(NA_character_)
  }
  normalizePath(existing[[1]], winslash = "/", mustWork = TRUE)
}

config_script_path <- find_first_existing(c(
  "utils/config.R",
  "../utils/config.R"
))

if (is.na(config_script_path)) {
  stop("Could not find utils/config.R from the current working directory.")
}

source(config_script_path)

config_path <- find_first_existing(c(
  "config/config.json",
  "../config/config.json"
))

if (is.na(config_path)) {
  stop("Could not find config/config.json from the current working directory.")
}

project_root <- dirname(dirname(config_path))
config <- load_config(config_path, required = TRUE)

out_dir <- get_config_value(
  config,
  "note_extraction_dir",
  default = file.path(project_root, "output", "note_extraction")
)

notes_dir <- file.path(project_root, "notes")
dir.create(notes_dir, recursive = TRUE, showWarnings = FALSE)

handp_file <- get_config_value(
  config,
  "handp_notes_file",
  default = file.path(notes_dir, "HP_NOTES.csv")
)

radiology_file <- get_config_value(
  config,
  "radiology_notes_file",
  default = file.path(notes_dir, "RAD_NOTES.csv")
)

dictionary_file <- file.path(project_root, "utils", "nlp_dictionary_concepts.csv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(dictionary_file)) {
  stop("Dictionary file not found: ", dictionary_file)
}

safe_read_csv_required <- function(path) {
  if (!file.exists(path)) {
    stop("Required input file not found: ", path)
  }
  readr::read_csv(path, show_col_types = FALSE)
}

normalize_text <- function(x) {
  x %>%
    replace_na("") %>%
    str_replace_all("[\r\n\t]+", " ") %>%
    str_squish() %>%
    str_to_lower()
}

parse_flexible_datetime <- function(x) {
  suppressWarnings(parse_date_time(
    x,
    orders = c(
      "ymd HMS", "ymd HM", "ymd",
      "mdy HMS", "mdy HM", "mdy",
      "Ymd HMS", "Ymd HM", "Ymd"
    ),
    tz = "America/Chicago"
  ))
}

column_or_na <- function(df, col_name) {
  if (col_name %in% names(df)) {
    return(df[[col_name]])
  }
  rep(NA_character_, nrow(df))
}

first_nonmissing <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) {
    return(NA_character_)
  }
  x[[1]]
}

extract_snippet <- function(text, pattern, window = 90) {
  if (is.na(text) || identical(text, "")) {
    return(NA_character_)
  }
  loc <- stringr::str_locate(text, regex(pattern, ignore_case = TRUE))
  if (all(is.na(loc))) {
    return(NA_character_)
  }
  start <- max(1, loc[1] - window)
  end <- min(nchar(text), loc[2] + window)
  str_sub(text, start, end) %>% str_squish()
}

collapse_patterns <- function(dict, domain, concept_group) {
  dict %>%
    filter(note_domain %in% c("all", domain), concept_group == !!concept_group) %>%
    pull(pattern) %>%
    unique() %>%
    paste(collapse = "|")
}

derive_rule_label <- function(df) {
  df %>%
    mutate(
      encounter_rule_label = case_when(
        encounter_any_ici_exposure == 1 &
          (encounter_any_irae_signal == 1 | encounter_radiology_supportive == 1) ~ "high_priority_ici_irae_review",
        encounter_any_irae_signal == 1 |
          (encounter_any_ici_exposure == 1 & (encounter_steroid_signal == 1 | encounter_rescue_signal == 1)) ~ "moderate_priority_possible_ici_irae",
        encounter_any_ici_exposure == 0 &
          encounter_any_irae_signal == 0 &
          (encounter_competing_infection == 1 | encounter_competing_progression == 1 | encounter_competing_edema == 1 | encounter_competing_aspiration == 1) ~ "competing_diagnosis_only",
        TRUE ~ "low_signal"
      ),
      encounter_llm_priority = as.integer(encounter_rule_label %in% c(
        "high_priority_ici_irae_review",
        "moderate_priority_possible_ici_irae"
      ))
    )
}

score_note <- function(df) {
  df %>%
    mutate(
      rule_score =
        3 * any_ici_agent +
        2 * any_irae_general +
        2 * any_irae_pneumonitis +
        2 * any_irae_myocarditis +
        2 * any_irae_hepatitis +
        2 * any_irae_colitis +
        2 * any_irae_nephritis +
        1 * any_irae_endocrine +
        1 * any_irae_neuro +
        1 * any_irae_hlh +
        1 * any_treatment_steroid +
        2 * any_treatment_rescue +
        2 * any_radiology_supportive_pneumonitis +
        1 * any_qualifier_favoring +
        0 * any_qualifier_uncertainty,
      llm_priority_note = as.integer(
        rule_score >= 3 |
          (candidate_ici_exposure_note == 1 & candidate_irae_note == 1) |
          (note_source == "radiology" & any_radiology_supportive_pneumonitis == 1)
      )
    )
}

ensure_dictionary_columns <- function(df, dict) {
  concept_groups <- dict %>% distinct(concept_group) %>% pull(concept_group)
  for (group in concept_groups) {
    any_col <- paste0("any_", group)
    snippet_col <- paste0("snippet_", group)
    if (!any_col %in% names(df)) {
      df[[any_col]] <- 0L
    }
    if (!snippet_col %in% names(df)) {
      df[[snippet_col]] <- NA_character_
    }
  }
  df
}

apply_dictionary_flags <- function(df, dict, domain) {
  concept_groups <- dict %>%
    filter(note_domain %in% c("all", domain)) %>%
    distinct(concept_group) %>%
    pull(concept_group)
  
  out <- df
  for (group in concept_groups) {
    patt <- collapse_patterns(dict, domain, group)
    col_name <- paste0("any_", group)
    out[[col_name]] <- as.integer(str_detect(out$note_text_clean, regex(patt, ignore_case = TRUE)))
    snippet_col <- paste0("snippet_", group)
    out[[snippet_col]] <- purrr::map_chr(out$NOTE_TEXT, extract_snippet, pattern = patt)
  }
  out
}

build_note_level <- function(df, dict, domain, note_source) {
  df %>%
    mutate(
      note_source = note_source,
      note_text_clean = normalize_text(NOTE_TEXT)
    ) %>%
    apply_dictionary_flags(dict, domain = domain) %>%
    ensure_dictionary_columns(dict = dict) %>%
    mutate(
      candidate_ici_exposure_note = as.integer(any_ici_agent == 1 | any_ici_general == 1),
      candidate_irae_note = as.integer(
        any_irae_general == 1 |
          any_irae_pneumonitis == 1 |
          any_irae_myocarditis == 1 |
          any_irae_hepatitis == 1 |
          any_irae_colitis == 1 |
          any_irae_nephritis == 1 |
          any_irae_endocrine == 1 |
          any_irae_neuro == 1 |
          any_irae_hlh == 1
      ),
      note_datetime = coalesce(
        parse_flexible_datetime(column_or_na(cur_data(), "NOTE_DTTM")),
        parse_flexible_datetime(column_or_na(cur_data(), "FINALIZING_DTTM")),
        parse_flexible_datetime(column_or_na(cur_data(), "SERVICE_DTTM"))
      )
    ) %>%
    score_note()
}

aggregate_encounter_level <- function(note_df) {
  note_df %>%
    group_by(MRN, HAR) %>%
    summarise(
      n_notes = n(),
      n_llm_priority_notes = sum(llm_priority_note, na.rm = TRUE),
      first_note_datetime = suppressWarnings(min(note_datetime, na.rm = TRUE)),
      last_note_datetime = suppressWarnings(max(note_datetime, na.rm = TRUE)),
      encounter_any_ici_exposure = as.integer(any(candidate_ici_exposure_note == 1, na.rm = TRUE)),
      encounter_explicit_ici_agent = as.integer(any(any_ici_agent == 1, na.rm = TRUE)),
      encounter_any_irae_signal = as.integer(any(candidate_irae_note == 1, na.rm = TRUE)),
      encounter_pneumonitis_signal = as.integer(any(any_irae_pneumonitis == 1, na.rm = TRUE)),
      encounter_cardiac_signal = as.integer(any(any_irae_myocarditis == 1, na.rm = TRUE)),
      encounter_hepatitis_signal = as.integer(any(any_irae_hepatitis == 1, na.rm = TRUE)),
      encounter_colitis_signal = as.integer(any(any_irae_colitis == 1, na.rm = TRUE)),
      encounter_renal_signal = as.integer(any(any_irae_nephritis == 1, na.rm = TRUE)),
      encounter_endocrine_signal = as.integer(any(any_irae_endocrine == 1, na.rm = TRUE)),
      encounter_neurologic_signal = as.integer(any(any_irae_neuro == 1, na.rm = TRUE)),
      encounter_hlh_signal = as.integer(any(any_irae_hlh == 1, na.rm = TRUE)),
      encounter_competing_infection = as.integer(any(any_competing_infection == 1, na.rm = TRUE)),
      encounter_competing_progression = as.integer(any(any_competing_progression == 1, na.rm = TRUE)),
      encounter_competing_edema = as.integer(any(any_competing_edema == 1, na.rm = TRUE)),
      encounter_competing_aspiration = as.integer(any(any_competing_aspiration == 1, na.rm = TRUE)),
      encounter_steroid_signal = as.integer(any(any_treatment_steroid == 1, na.rm = TRUE)),
      encounter_rescue_signal = as.integer(any(any_treatment_rescue == 1, na.rm = TRUE)),
      encounter_radiology_supportive = as.integer(any(any_radiology_supportive_pneumonitis == 1, na.rm = TRUE)),
      max_rule_score = suppressWarnings(max(rule_score, na.rm = TRUE)),
      example_ici_snippet = first_nonmissing(c(snippet_ici_agent, snippet_ici_general)),
      example_irae_snippet = first_nonmissing(c(
        snippet_irae_general, snippet_irae_pneumonitis, snippet_irae_myocarditis,
        snippet_irae_hepatitis, snippet_irae_colitis, snippet_irae_nephritis,
        snippet_irae_endocrine, snippet_irae_neuro, snippet_irae_hlh
      )),
      example_competing_snippet = first_nonmissing(c(
        snippet_competing_infection, snippet_competing_progression,
        snippet_competing_edema, snippet_competing_aspiration
      )),
      .groups = "drop"
    ) %>%
    mutate(
      first_note_datetime = if_else(is.infinite(first_note_datetime), as.POSIXct(NA), first_note_datetime),
      last_note_datetime = if_else(is.infinite(last_note_datetime), as.POSIXct(NA), last_note_datetime)
    ) %>%
    derive_rule_label()
}

prepare_llm_queue <- function(note_df, encounter_df, max_handp = 2, max_radiology = 3) {
  encounter_priority <- encounter_df %>%
    filter(encounter_llm_priority == 1) %>%
    select(MRN, HAR, encounter_rule_label)
  
  note_df %>%
    inner_join(encounter_priority, by = c("MRN", "HAR")) %>%
    arrange(MRN, HAR, desc(rule_score), note_datetime) %>%
    group_by(MRN, HAR, note_source) %>%
    mutate(note_rank_within_source = row_number()) %>%
    ungroup() %>%
    filter(
      (note_source == "handp" & note_rank_within_source <= max_handp) |
        (note_source == "radiology" & note_rank_within_source <= max_radiology)
    ) %>%
    select(any_of(c(
      MRN, HAR, note_source, NOTE_ID, NOTE_CSN_ID, note_datetime,
      NOTE_TYPE, PROC_NAME, AUTHOR_SERV, rule_score, llm_priority_note,
      encounter_rule_label,
      snippet_ici_agent, snippet_ici_general,
      snippet_irae_general, snippet_irae_pneumonitis, snippet_irae_myocarditis,
      snippet_competing_infection, snippet_competing_progression,
      snippet_radiology_supportive_pneumonitis,
      NOTE_TEXT
    )))
}

dictionary <- safe_read_csv_required(dictionary_file)

message("Looking for H&P notes at: ", handp_file)
message("Looking for radiology reports at: ", radiology_file)

if (!file.exists(handp_file) && !file.exists(radiology_file)) {
  stop(
    "No note source files found. Place `HP_NOTES.csv` and/or `RAD_NOTES.csv` ",
    "in the repo `notes/` directory, or specify `handp_notes_file` and/or ",
    "`radiology_notes_file` in config/config.json."
  )
}

note_tables <- list()

if (file.exists(handp_file)) {
  handp_raw <- safe_read_csv_required(handp_file)
  required_handp <- c("MRN", "HAR", "NOTE_ID", "NOTE_TEXT")
  missing_handp <- setdiff(required_handp, names(handp_raw))
  if (length(missing_handp) > 0) {
    stop("H&P file missing required columns: ", paste(missing_handp, collapse = ", "))
  }
  
  handp_note_level <- build_note_level(handp_raw, dictionary, domain = "all", note_source = "handp")
  write_csv(handp_note_level, file.path(out_dir, "handp_note_level_rules.csv"))
  note_tables$handp <- handp_note_level
}

if (file.exists(radiology_file)) {
  radiology_raw <- safe_read_csv_required(radiology_file)
  required_img <- c("MRN", "HAR", "NOTE_ID", "NOTE_TEXT", "PROC_NAME")
  missing_img <- setdiff(required_img, names(radiology_raw))
  if (length(missing_img) > 0) {
    stop("Radiology file missing required columns: ", paste(missing_img, collapse = ", "))
  }
  
  chest_proc_pattern <- "chest|lung|thorax|pe chest|cta pe|ild"
  radiology_filtered <- radiology_raw %>%
    mutate(proc_name_clean = normalize_text(PROC_NAME)) %>%
    filter(str_detect(proc_name_clean, chest_proc_pattern) | str_detect(normalize_text(NOTE_TEXT), "chest|lung|pulmonary|pneumon"))
  
  radiology_note_level <- build_note_level(radiology_filtered, dictionary, domain = "radiology", note_source = "radiology")
  write_csv(radiology_note_level, file.path(out_dir, "radiology_note_level_rules.csv"))
  note_tables$radiology <- radiology_note_level
}

all_note_level <- bind_rows(note_tables) %>%
  ensure_dictionary_columns(dict = dictionary)
write_csv(all_note_level, file.path(out_dir, "all_note_level_rules.csv"))

encounter_level <- aggregate_encounter_level(all_note_level)
write_csv(encounter_level, file.path(out_dir, "encounter_level_rule_labels.csv"))

llm_queue <- prepare_llm_queue(all_note_level, encounter_level)
write_csv(llm_queue, file.path(out_dir, "llm_review_queue.csv"))

summary_tbl <- tibble(
  n_notes_total = nrow(all_note_level),
  n_handp_notes = sum(all_note_level$note_source == "handp"),
  n_radiology_notes = sum(all_note_level$note_source == "radiology"),
  n_encounters = nrow(encounter_level),
  n_high_priority = sum(encounter_level$encounter_rule_label == "high_priority_ici_irae_review"),
  n_moderate_priority = sum(encounter_level$encounter_rule_label == "moderate_priority_possible_ici_irae"),
  n_llm_queue_notes = nrow(llm_queue)
)

write_csv(summary_tbl, file.path(out_dir, "extraction_summary.csv"))
print(summary_tbl)
