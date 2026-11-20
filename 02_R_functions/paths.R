# paths.R
# Path helpers on top of here::here(). The names of 03_input_files, 04_output,
# 05_manuscript and 08_archive are defined here, for every path that a step
# builds through these helpers. Other places also name a numbered folder and
# must be edited with this file if one is renamed: the line of every script that
# sources 02_R_functions/, run_pipeline.R (01_pipeline), the two build scripts
# under 03_input_files/, step 13 and the scripts in 06_diagnostics/ (04_output
# and 06_diagnostics), the source column of
# 03_input_files/manuscript/manuscript_items_map.csv, two expressions of
# manuscript_values_map.csv in the same folder (rows D4.8-002 and R3.4-007 name
# 03_input_files), and the patterns of .gitignore and .gitattributes. The root
# README lists them under "How paths work".
#
# Write guard. In an optional run (config$run_kind "optional": an optional
# analysis, the gene-level GO enrichment or the diagnostics turned on) nothing
# may be written into a folder of the paper. output_dir() and output_path()
# therefore stop in such a run, except for a step's cache folder and for folders
# under 04_output/not_adopted/; an optional analysis writes through
# analysis_dir() and analysis_path() (02_R_functions/analyses.R), and reads the
# paper's outputs through paper_path(). manuscript_path() allows only
# 05_manuscript/not_adopted/ in such a run, and never a ".." component.

input_path <- function(...) here::here("03_input_files", ...)

optional_run <- function() exists("config", inherits = TRUE) && identical(config$run_kind, "optional")

# Folder for one step's outputs, created if absent: 04_output/<step>/<sub>/
output_dir <- function(step, sub = NULL) {
  if (optional_run() && !identical(sub, "cache") && !startsWith(step, "not_adopted")) {
    stop("An optional run may not write into the paper folder 04_output/", step, if (!is.null(sub)) paste0("/", sub),
         "; use analysis_path().", call. = FALSE)
  }
  d <- if (is.null(sub)) here::here("04_output", step) else here::here("04_output", step, sub)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

output_path <- function(step, sub, file) file.path(output_dir(step, sub), file)

# A file of the paper's outputs, for reading (never creates a folder).
paper_path <- function(step, sub, file) here::here("04_output", step, sub, file)

# Folder for a step's large, regenerable intermediates (git-ignored).
step_cache_dir <- function(step) output_dir(step, "cache")

manuscript_path <- function(...) {
  parts <- as.character(c(...))
  if (any(unlist(strsplit(parts, "[/\\\\]")) == "..")) stop("A path under 05_manuscript/ may not contain '..': ", paste(parts, collapse = "/"), call. = FALSE)
  if (optional_run() && (length(parts) == 0 || parts[1] != "not_adopted")) {
    stop("An optional run may write into 05_manuscript/not_adopted/ only.", call. = FALSE)
  }
  here::here("05_manuscript", ...)
}
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
