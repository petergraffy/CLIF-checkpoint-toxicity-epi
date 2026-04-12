suppressPackageStartupMessages({
  library(jsonlite)
})

load_config <- function(json_path = "config/config.json", required = TRUE) {
  if (!file.exists(json_path)) {
    if (required) {
      stop(
        "Configuration file not found at ", json_path, ". ",
        "Create config/config.json from config/config_template.json."
      )
    }
    return(list())
  }
  
  config <- jsonlite::fromJSON(json_path, simplifyVector = TRUE)
  if (is.null(config$clif_dir) && !is.null(config$tables_path)) {
    config$clif_dir <- config$tables_path
  }
  message("Loaded configuration from ", json_path)
  config
}

get_config_value <- function(config, key, default = NULL, required = FALSE) {
  value <- config[[key]]
  if (is.null(value) || identical(value, "")) {
    if (required) {
      stop("Missing required config field: ", key)
    }
    return(default)
  }
  value
}

config <- tryCatch(
  load_config(required = FALSE),
  error = function(e) {
    message("Config not loaded: ", e$message)
    list()
  }
)
