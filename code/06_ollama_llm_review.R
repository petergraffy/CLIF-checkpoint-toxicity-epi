# ============================================================
# 06_ollama_llm_review.R
# Run local Ollama review over H&P encounter packets
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
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

encounter_packets_path <- file.path(llm_review_dir, "handp_llm_encounter_packets.csv")
manual_review_template_path <- file.path(llm_review_dir, "handp_manual_review_template.csv")

ollama_model <- get_config_value(config, "ollama_model", default = "qwen3:235b")
ollama_host <- get_config_value(config, "ollama_host", default = "http://127.0.0.1:11434")
ollama_cli <- get_config_value(config, "ollama_cli", default = Sys.which("ollama"))
ollama_mode <- get_config_value(config, "ollama_mode", default = "auto")
ollama_temperature <- as.numeric(get_config_value(config, "ollama_temperature", default = 0))
ollama_max_packets <- as.integer(get_config_value(config, "ollama_max_packets", default = NA))
ollama_sleep_seconds <- as.numeric(get_config_value(config, "ollama_sleep_seconds", default = 0))

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

normalize_json_text <- function(x) {
  x <- sanitize_character_encoding(x)
  x <- trimws(x)
  x <- sub("^```json\\s*", "", x, perl = TRUE)
  x <- sub("^```\\s*", "", x, perl = TRUE)
  x <- sub("\\s*```$", "", x, perl = TRUE)
  trimws(x)
}

detect_ollama_runtime <- function(mode = "auto", cli_path = "", host = "http://127.0.0.1:11434") {
  cli_available <- nzchar(cli_path) && file.exists(cli_path)
  api_available <- FALSE

  if (mode %in% c("auto", "api")) {
    probe <- suppressWarnings(
      system2(
        "curl",
        c("-s", paste0(host, "/api/tags")),
        stdout = TRUE,
        stderr = TRUE
      )
    )
    api_available <- length(probe) > 0 && any(grepl("models", probe, fixed = TRUE))
  }

  if (mode == "cli" && cli_available) {
    return("cli")
  }
  if (mode == "api" && api_available) {
    return("api")
  }
  if (mode == "auto" && cli_available) {
    return("cli")
  }
  if (mode == "auto" && api_available) {
    return("api")
  }

  NA_character_
}

call_ollama_cli <- function(prompt, model, cli_path) {
  out <- system2(
    cli_path,
    c("run", model),
    input = prompt,
    stdout = TRUE,
    stderr = TRUE
  )
  paste(out, collapse = "\n")
}

call_ollama_api <- function(prompt, model, host, temperature = 0) {
  payload <- list(
    model = model,
    prompt = prompt,
    stream = FALSE,
    options = list(
      temperature = temperature
    )
  )

  payload_path <- tempfile(fileext = ".json")
  writeLines(
    jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null"),
    con = payload_path,
    useBytes = TRUE
  )

  out <- system2(
    "curl",
    c(
      "-s",
      "-X", "POST",
      "-H", "Content-Type: application/json",
      paste0(host, "/api/generate"),
      "-d", paste0("@", payload_path)
    ),
    stdout = TRUE,
    stderr = TRUE
  )

  unlink(payload_path)
  paste(out, collapse = "\n")
}

extract_response_text <- function(runtime_mode, raw_response) {
  if (runtime_mode == "cli") {
    return(normalize_json_text(raw_response))
  }

  parsed <- tryCatch(
    jsonlite::fromJSON(raw_response, simplifyVector = FALSE),
    error = function(e) NULL
  )

  if (is.null(parsed) || is.null(parsed$response)) {
    return(normalize_json_text(raw_response))
  }

  normalize_json_text(parsed$response)
}

parse_llm_json <- function(text) {
  parsed <- tryCatch(
    jsonlite::fromJSON(text, simplifyVector = FALSE),
    error = function(e) NULL
  )

  if (is.null(parsed)) {
    return(NULL)
  }

  parsed
}

collapse_json_field <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(NA_character_)
  }
  if (is.list(x) && !is.atomic(x)) {
    x <- unlist(x, use.names = FALSE)
  }
  vals <- as.character(x)
  vals <- vals[!is.na(vals) & vals != ""]
  if (length(vals) == 0) {
    return(NA_character_)
  }
  paste(unique(vals), collapse = " | ")
}

flatten_llm_record <- function(packet_row, parsed_json, response_text, parse_ok, error_message = NA_character_) {
  tibble(
    TEMP_ID = packet_row$TEMP_ID,
    MRN = packet_row$MRN,
    HAR = packet_row$HAR,
    encounter_rule_label = packet_row$encounter_rule_label,
    llm_model = ollama_model,
    parse_ok = as.integer(parse_ok),
    llm_error = error_message,
    llm_response_text = response_text,
    ici_exposure_status = if (parse_ok) collapse_json_field(parsed_json$ici_exposure_status) else NA_character_,
    ici_agents_mentioned = if (parse_ok) collapse_json_field(parsed_json$ici_agents_mentioned) else NA_character_,
    ici_recent_or_active = if (parse_ok) collapse_json_field(parsed_json$ici_recent_or_active) else NA_character_,
    suspected_irae_status = if (parse_ok) collapse_json_field(parsed_json$suspected_irae_status) else NA_character_,
    irae_organ_systems = if (parse_ok) collapse_json_field(parsed_json$irae_organ_systems) else NA_character_,
    confidence_irae_primary_driver = if (parse_ok) collapse_json_field(parsed_json$confidence_irAE_is_primary_driver) else NA_character_,
    competing_diagnoses = if (parse_ok) collapse_json_field(parsed_json$competing_diagnoses) else NA_character_,
    steroids_for_suspected_irae = if (parse_ok) collapse_json_field(parsed_json$steroids_for_suspected_irae) else NA_character_,
    rescue_immunosuppression_for_suspected_irae = if (parse_ok) collapse_json_field(parsed_json$rescue_immunosuppression_for_suspected_irae) else NA_character_,
    overall_llm_label = if (parse_ok) collapse_json_field(parsed_json$overall_llm_label) else NA_character_,
    supporting_evidence = if (parse_ok) collapse_json_field(parsed_json$supporting_evidence) else NA_character_,
    reasoning_summary = if (parse_ok) collapse_json_field(parsed_json$reasoning_summary) else NA_character_
  )
}

if (!file.exists(encounter_packets_path)) {
  stop(
    "Could not find handp_llm_encounter_packets.csv in ", llm_review_dir, ". ",
    "Run code/05_llm_review_prep.R first."
  )
}

runtime_mode <- detect_ollama_runtime(
  mode = ollama_mode,
  cli_path = ollama_cli,
  host = ollama_host
)

if (is.na(runtime_mode)) {
  stop(
    "No Ollama runtime detected. Install the `ollama` CLI or start a local API at ",
    ollama_host, "."
  )
}

encounter_packets <- safe_read_csv_required(encounter_packets_path) %>%
  mutate(
    HAR = as.character(HAR),
    MRN = sanitize_character_encoding(MRN),
    prompt_text = sanitize_character_encoding(prompt_text),
    encounter_rule_label = sanitize_character_encoding(encounter_rule_label)
  )

if (!is.na(ollama_max_packets)) {
  encounter_packets <- encounter_packets %>%
    slice_head(n = ollama_max_packets)
}

if (nrow(encounter_packets) == 0) {
  stop("No encounter packets available for Ollama review.")
}

raw_results <- vector("list", nrow(encounter_packets))
flat_results <- vector("list", nrow(encounter_packets))

for (i in seq_len(nrow(encounter_packets))) {
  packet_row <- encounter_packets[i, ]
  prompt <- packet_row$prompt_text[[1]]

  raw_response <- tryCatch(
    {
      if (runtime_mode == "cli") {
        call_ollama_cli(prompt, model = ollama_model, cli_path = ollama_cli)
      } else {
        call_ollama_api(prompt, model = ollama_model, host = ollama_host, temperature = ollama_temperature)
      }
    },
    error = function(e) structure(e$message, class = "ollama_error")
  )

  if (inherits(raw_response, "ollama_error")) {
    raw_text <- as.character(raw_response)
    flat_results[[i]] <- flatten_llm_record(
      packet_row = packet_row,
      parsed_json = NULL,
      response_text = raw_text,
      parse_ok = FALSE,
      error_message = raw_text
    )
    raw_results[[i]] <- tibble(
      TEMP_ID = packet_row$TEMP_ID,
      HAR = packet_row$HAR,
      llm_model = ollama_model,
      runtime_mode = runtime_mode,
      raw_response = raw_text
    )
    next
  }

  response_text <- extract_response_text(runtime_mode, raw_response)
  parsed_json <- parse_llm_json(response_text)
  parse_ok <- !is.null(parsed_json)

  flat_results[[i]] <- flatten_llm_record(
    packet_row = packet_row,
    parsed_json = parsed_json,
    response_text = response_text,
    parse_ok = parse_ok,
    error_message = if (parse_ok) NA_character_ else "Failed to parse LLM response as JSON"
  )

  raw_results[[i]] <- tibble(
    TEMP_ID = packet_row$TEMP_ID,
    HAR = packet_row$HAR,
    llm_model = ollama_model,
    runtime_mode = runtime_mode,
    raw_response = response_text
  )

  if (ollama_sleep_seconds > 0) {
    Sys.sleep(ollama_sleep_seconds)
  }
}

raw_results_tbl <- bind_rows(raw_results)
flat_results_tbl <- bind_rows(flat_results)

review_merged <- if (file.exists(manual_review_template_path)) {
  manual_review <- safe_read_csv_required(manual_review_template_path) %>%
    mutate(HAR = as.character(HAR))

  manual_review %>%
    left_join(
      flat_results_tbl %>%
        select(
          HAR,
          llm_model,
          parse_ok,
          llm_error,
          ici_exposure_status,
          ici_agents_mentioned,
          ici_recent_or_active,
          suspected_irae_status,
          irae_organ_systems,
          confidence_irae_primary_driver,
          competing_diagnoses,
          steroids_for_suspected_irae,
          rescue_immunosuppression_for_suspected_irae,
          overall_llm_label,
          supporting_evidence,
          reasoning_summary
        ),
      by = "HAR"
    )
} else {
  flat_results_tbl
}

summary_tbl <- tibble(
  runtime_mode = runtime_mode,
  llm_model = ollama_model,
  n_packets_attempted = nrow(encounter_packets),
  n_parse_ok = sum(flat_results_tbl$parse_ok, na.rm = TRUE),
  n_parse_failed = sum(flat_results_tbl$parse_ok == 0, na.rm = TRUE),
  n_likely_ici_irae = sum(flat_results_tbl$overall_llm_label == "likely_ici_irae", na.rm = TRUE),
  n_possible_ici_irae = sum(flat_results_tbl$overall_llm_label == "possible_ici_irae", na.rm = TRUE),
  n_unlikely_ici_irae = sum(flat_results_tbl$overall_llm_label == "unlikely_ici_irae", na.rm = TRUE)
)

write_csv(raw_results_tbl, file.path(llm_review_dir, "ollama_handp_raw_responses.csv"))
write_csv(flat_results_tbl, file.path(llm_review_dir, "ollama_handp_parsed_outputs.csv"))
write_csv(review_merged, file.path(llm_review_dir, "ollama_handp_review_merged.csv"))
write_csv(summary_tbl, file.path(llm_review_dir, "ollama_handp_review_summary.csv"))

print(summary_tbl)
