# Shared input/output helper functions for the gene-X SCI bulk RNA-seq analysis.
# These functions keep path handling consistent across all R scripts.

suppressPackageStartupMessages({
  library(yaml)
})

load_paths <- function(config_file = "config/paths.yml") {
  if (!file.exists(config_file)) {
    stop("Missing config file: ", config_file,
         "\nCopy config/paths_template.yml to config/paths.yml and edit local paths.")
  }
  yaml::read_yaml(config_file)
}

project_path <- function(..., paths = NULL) {
  if (is.null(paths)) {
    paths <- load_paths()
  }
  file.path(paths$project$root, ...)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

check_file_exists <- function(path, label = NULL) {
  if (!file.exists(path)) {
    if (is.null(label)) label <- path
    stop("Missing required file: ", label, "\nExpected path: ", path)
  }
  invisible(path)
}
