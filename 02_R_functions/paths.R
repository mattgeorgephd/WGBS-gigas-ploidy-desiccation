# paths.R
# Path helpers on top of here::here(). The names of 03_input_files, 04_output,
# 05_manuscript and 08_archive are defined here, for every path that a step
# builds through these helpers. Other places also name a numbered folder and
# must be edited with this file if one is renamed: the line of every script that
# sources 02_R_functions/, run_pipeline.R (01_pipeline), the two build scripts
# under 03_input_files/, step 13 and the scripts in 06_diagnostics/ (04_output
# and 06_diagnostics), and the source column of
# 03_input_files/manuscript/manuscript_items_map.csv. The root README lists them
# under "How paths work".

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
# No step of the pipeline reads the archive: the legacy files that steps 11 and 12
# compare with are copies under 03_input_files/. archive_path() serves the checks
# in 06_diagnostics/, which compare pipeline outputs with the legacy outputs.
archive_path    <- function(...) here::here("08_archive", ...)

# Folder that holds the files downloaded from gannet. `config$cache_dir` (set
# from the environment variable WGBS_CACHE_DIR) overrides the default.
gannet_cache_dir <- function(config, sub = NULL) {
  base <- if (nzchar(config$cache_dir)) config$cache_dir else input_path("gannet_cache")
  d <- if (is.null(sub)) base else file.path(base, sub)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
