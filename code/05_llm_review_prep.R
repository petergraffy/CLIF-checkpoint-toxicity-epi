# ============================================================
# 05_llm_review_prep.R
# Prepare H&P-focused encounter packets for LLM review/adjudication
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(jsonlite)
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

note_extraction_dir <- get_config_value(
  config,
  "note_extraction_dir",
  default = file.path(project_root, "output", "note_extraction")
)

llm_review_dir <- file.path(note_extraction_dir, "llm_review")
dir.create(llm_review_dir, recursive = TRUE, showWarnings = FALSE)

analysis_dir <- get_config_value(
  config,
  "analysis_dir",
  default = file.path(project_root, "output", "checkpoint_irae_icu_epi")
)

analysis_path <- file.path(analysis_dir, "analysis_dataset.csv")
encounter_labels_path <- file.path(note_extraction_dir, "encounter_level_rule_labels.csv")
all_priority_notes_path <- file.path(note_extraction_dir, "all_notes_high_moderate_encounters.csv")
strict_priority_notes_path <- file.path(note_extraction_dir, "high_moderate_priority_notes_only.csv")

safe_read_csv_required <- function(path) {
  if (!file.exists(path)) {
    stop("Required input file not found: ", path)
  }

  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  raw_head <- readBin(con, what = "raw", n = 4000)

  is_utf16le_bom <- length(raw_head) >= 2 && identical(as.integer(raw_head[1:2]), c(255L, 254L))
  is_utf16be_bom <- length(raw_head) >= 2 && identical(as.integer(raw_head[1:2]), c(254L, 255L))
  has_embedded_nul <- any(raw_head == as.raw(0))

  if (is_utf16le_bom || has_embedded_nul) {
    return(readr::read_csv(path, locale = readr::locale(encoding = "UTF-16LE"), show_col_types = FALSE))
  }

  if (is_utf16be_bom) {
    return(readr::read_csv(path, locale = readr::locale(encoding = "UTF-16BE"), show_col_types = FALSE))
  }

  readr::read_csv(path, show_col_types = FALSE)
}

sanitize_character_encoding <- function(x) {
  x <- as.character(x)
  out <- iconv(x, from = "", to = "UTF-8", sub = " ")
  fallback <- is.na(out) & !is.na(x)

  if (any(fallback)) {
    out[fallback] <- tryCatch(
      iconv(x[fallback], from = "CP1252", to = "UTF-8", sub = " "),
      error = function(e) iconv(x[fallback], from = "latin1", to = "UTF-8", sub = " ")
    )
  }

  out
}

normalize_text <- function(x) {
  sanitize_character_encoding(x) %>%
    replace_na("") %>%
    str_replace_all("[\r\n\t]+", " ") %>%
    str_squish()
}

first_nonmissing <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) {
    return(NA_character_)
  }
  x[[1]]
}

collapse_unique <- function(x, sep = " | ") {
  vals <- unique(x[!is.na(x) & x != "" & x != "NA"])
  if (length(vals) == 0) {
    return(NA_character_)
  }
  paste(vals, collapse = sep)
}

to_binary_flag <- function(x) {
  as.integer(!is.na(x) & x != "" & x != "NA")
}

build_note_excerpt <- function(note_row, note_index, max_chars = 2200) {
  parts <- c(
    paste0("### Note ", note_index),
    paste0("Note source: ", note_row[["note_source"]]),
    paste0("Note datetime: ", note_row[["note_datetime"]]),
    paste0("Note type: ", note_row[["NOTE_TYPE"]]),
    paste0("Author service: ", note_row[["AUTHOR_SERV"]]),
    paste0("Rule score: ", note_row[["rule_score"]]),
    "Note text:",
    str_sub(note_row[["NOTE_TEXT"]], 1, max_chars)
  )
  paste(parts, collapse = "\n")
}

build_encounter_packet <- function(note_df) {
  paste(
    purrr::map_chr(seq_len(nrow(note_df)), function(idx) {
      build_note_excerpt(as.list(note_df[idx, ]), note_index = idx)
    }),
    collapse = "\n\n"
  )
}

build_llm_prompt <- function(encounter_row) {
  paste(
    "You are reviewing hospital H&P notes for possible immune checkpoint inhibitor toxicity.",
    "Use only the provided note text.",
    "Return strict JSON with the schema below and do not add commentary outside JSON.",
    "",
    "Required JSON fields:",
    "{",
    '  "ici_exposure_status": "yes|no|uncertain",',
    '  "ici_agents_mentioned": ["string"],',
    '  "ici_recent_or_active": "yes|no|uncertain",',
    '  "suspected_irae_status": "yes|no|uncertain",',
    '  "irae_organ_systems": ["pneumonitis|myocarditis|hepatitis|colitis|nephritis|endocrine|neurologic|hlh_like|other"],',
    '  "confidence_irAE_is_primary_driver": "high|moderate|low",',
    '  "competing_diagnoses": ["infection|progression|edema|aspiration|other"],',
    '  "steroids_for_suspected_irae": "yes|no|uncertain",',
    '  "rescue_immunosuppression_for_suspected_irae": "yes|no|uncertain",',
    '  "overall_llm_label": "likely_ici_irae|possible_ici_irae|unlikely_ici_irae",',
    '  "supporting_evidence": ["short quoted snippets"],',
    '  "reasoning_summary": "1-3 sentence summary"',
    "}",
    "",
    "Encounter metadata:",
    paste0("TEMP_ID: ", encounter_row[["TEMP_ID"]]),
    paste0("Encounter rule label: ", encounter_row[["encounter_rule_label"]]),
    paste0("Top note count: ", encounter_row[["n_notes_in_packet"]]),
    "",
    "Notes:",
    encounter_row[["packet_text"]],
    sep = "\n"
  )
}

if (!file.exists(all_priority_notes_path)) {
  stop(
    "Could not find all_notes_high_moderate_encounters.csv in ", note_extraction_dir, ". ",
    "Run code/04_note_extraction.R first."
  )
}

if (!file.exists(encounter_labels_path)) {
  stop(
    "Could not find encounter_level_rule_labels.csv in ", note_extraction_dir, ". ",
    "Run code/04_note_extraction.R first."
  )
}

encounter_labels <- safe_read_csv_required(encounter_labels_path) %>%
  mutate(HAR = as.character(HAR))

analysis_subset <- if (file.exists(analysis_path)) {
  safe_read_csv_required(analysis_path) %>%
    transmute(
      HAR = as.character(hospitalization_id),
      patient_id = as.character(patient_id),
      hospitalization_id = as.character(hospitalization_id),
      first_icu_in,
      discharge_dttm,
      icu_year,
      phenotype_any,
      phenotype_primary,
      phenotype_high,
      probable_checkpoint_irae_icu,
      possible_checkpoint_irae_icu,
      high_confidence_checkpoint_irae_icu,
      cancer_type_broad
    ) %>%
    distinct(HAR, .keep_all = TRUE)
} else {
  tibble(HAR = character())
}

candidate_notes <- safe_read_csv_required(all_priority_notes_path) %>%
  mutate(
    HAR = as.character(HAR),
    note_source = sanitize_character_encoding(note_source),
    NOTE_TYPE = sanitize_character_encoding(NOTE_TYPE),
    AUTHOR_SERV = sanitize_character_encoding(AUTHOR_SERV),
    NOTE_TEXT = sanitize_character_encoding(NOTE_TEXT),
    rule_score = suppressWarnings(as.numeric(rule_score)),
    llm_priority_note = suppressWarnings(as.integer(llm_priority_note)),
    encounter_rule_label = sanitize_character_encoding(encounter_rule_label)
  ) %>%
  filter(note_source == "handp") %>%
  mutate(
    ici_snippet_present = to_binary_flag(snippet_ici_agent) | to_binary_flag(snippet_ici_general),
    irae_snippet_present = to_binary_flag(snippet_irae_general) | to_binary_flag(snippet_irae_pneumonitis) | to_binary_flag(snippet_irae_myocarditis),
    competing_snippet_present = to_binary_flag(snippet_competing_infection) | to_binary_flag(snippet_competing_progression),
    note_priority_rank = case_when(
      llm_priority_note == 1 ~ 1L,
      irae_snippet_present == 1 ~ 2L,
      ici_snippet_present == 1 ~ 3L,
      competing_snippet_present == 1 ~ 4L,
      TRUE ~ 5L
    )
  )

if (nrow(candidate_notes) == 0) {
  stop("No H&P notes were available in all_notes_high_moderate_encounters.csv.")
}

top_notes_per_encounter <- candidate_notes %>%
  arrange(HAR, note_priority_rank, desc(rule_score), note_datetime) %>%
  group_by(HAR) %>%
  mutate(note_rank_within_encounter = row_number()) %>%
  filter(note_rank_within_encounter <= 3) %>%
  ungroup()

encounter_packets <- top_notes_per_encounter %>%
  group_by(HAR) %>%
  summarise(
    encounter_rule_label = first_nonmissing(encounter_rule_label),
    n_notes_in_packet = n(),
    max_rule_score = suppressWarnings(max(rule_score, na.rm = TRUE)),
    any_note_llm_priority = as.integer(any(llm_priority_note == 1, na.rm = TRUE)),
    any_ici_snippet = as.integer(any(ici_snippet_present == 1, na.rm = TRUE)),
    any_irae_snippet = as.integer(any(irae_snippet_present == 1, na.rm = TRUE)),
    any_competing_snippet = as.integer(any(competing_snippet_present == 1, na.rm = TRUE)),
    evidence_ici = collapse_unique(c(snippet_ici_agent, snippet_ici_general)),
    evidence_irae = collapse_unique(c(snippet_irae_general, snippet_irae_pneumonitis, snippet_irae_myocarditis)),
    evidence_competing = collapse_unique(c(snippet_competing_infection, snippet_competing_progression)),
    packet_text = build_encounter_packet(cur_data()),
    .groups = "drop"
  ) %>%
  mutate(TEMP_ID = row_number()) %>%
  left_join(
    encounter_labels %>%
      select(
        HAR,
        MRN,
        n_notes,
        n_llm_priority_notes,
        encounter_any_ici_exposure,
        encounter_explicit_ici_agent,
        encounter_any_irae_signal,
        encounter_pneumonitis_signal,
        encounter_cardiac_signal,
        encounter_hepatitis_signal,
        encounter_colitis_signal,
        encounter_renal_signal,
        encounter_endocrine_signal,
        encounter_neurologic_signal,
        encounter_hlh_signal,
        encounter_competing_infection,
        encounter_competing_progression,
        encounter_steroid_signal,
        encounter_rescue_signal
      ),
    by = "HAR"
  ) %>%
  left_join(analysis_subset, by = "HAR") %>%
  rowwise() %>%
  mutate(prompt_text = build_llm_prompt(as.list(cur_data()))) %>%
  ungroup() %>%
  arrange(desc(any_note_llm_priority), desc(any_irae_snippet), desc(max_rule_score), HAR)

manual_review_template <- encounter_packets %>%
  transmute(
    TEMP_ID,
    MRN,
    HAR,
    patient_id,
    hospitalization_id,
    encounter_rule_label,
    icu_year,
    cancer_type_broad,
    n_notes_in_packet,
    max_rule_score,
    evidence_ici,
    evidence_irae,
    evidence_competing,
    reviewer_ici_exposure_status = NA_character_,
    reviewer_ici_agents = NA_character_,
    reviewer_ici_recent_or_active = NA_character_,
    reviewer_suspected_irae_status = NA_character_,
    reviewer_irae_organ_systems = NA_character_,
    reviewer_primary_competing_diagnosis = NA_character_,
    reviewer_steroids_for_suspected_irae = NA_character_,
    reviewer_rescue_immunosuppression = NA_character_,
    reviewer_overall_label = NA_character_,
    reviewer_confidence = NA_character_,
    reviewer_notes = NA_character_
  )

llm_jsonl <- encounter_packets %>%
  transmute(
    custom_id = paste0("har_", HAR),
    input = prompt_text
  )

write_csv(top_notes_per_encounter, file.path(llm_review_dir, "handp_top_notes_for_llm_review.csv"))
write_csv(encounter_packets, file.path(llm_review_dir, "handp_llm_encounter_packets.csv"))
write_csv(manual_review_template, file.path(llm_review_dir, "handp_manual_review_template.csv"))
write_lines(
  purrr::map_chr(seq_len(nrow(llm_jsonl)), function(i) {
    jsonlite::toJSON(llm_jsonl[i, ], auto_unbox = TRUE, null = "null")
  }),
  file.path(llm_review_dir, "handp_llm_prompts.jsonl")
)

summary_tbl <- tibble(
  n_handp_candidate_notes = nrow(candidate_notes),
  n_handp_top_notes = nrow(top_notes_per_encounter),
  n_handp_encounter_packets = nrow(encounter_packets),
  n_encounters_with_any_ici_snippet = sum(encounter_packets$any_ici_snippet, na.rm = TRUE),
  n_encounters_with_any_irae_snippet = sum(encounter_packets$any_irae_snippet, na.rm = TRUE),
  n_encounters_with_any_note_llm_priority = sum(encounter_packets$any_note_llm_priority, na.rm = TRUE)
)

write_csv(summary_tbl, file.path(llm_review_dir, "llm_review_prep_summary.csv"))
print(summary_tbl)
