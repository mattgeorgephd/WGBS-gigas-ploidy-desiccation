# paths.R
# Path helpers. Every path in the pipeline is built here, on top of here::here(),
# so the numbered folder names appear in exactly one file.

input_path <- function(...) here::here("03_input_files", ...)

# Folder for one step's outputs, created if absent: 04_output/<step>/<sub>/
output_dir <- function(step, sub = NULL) {
  d <- if (is.null(sub)) here::here("04_output", step) else here::here("04_output", step, sub)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

output_path <- function(step, sub, file) file.path(output_dir(step, sub), file)

# Folder for a step's large, regenerable intermediates (git-ignored).
step_cache_dir <- function(step) output_dir(step, "cache")

manuscript_path <- function(...) here::here("05_manuscript", ...)
archive_path    <- function(...) here::here("08_archive", ...)

# Folder that holds the files downloaded from gannet. `config$cache_dir` (set
# from the environment variable WGBS_CACHE_DIR) overrides the default.
gannet_cache_dir <- function(config, sub = NULL) {
  base <- if (nzchar(config$cache_dir)) config$cache_dir else input_path("gannet_cache")
  d <- if (is.null(sub)) base else file.path(base, sub)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
