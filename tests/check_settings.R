# tests/check_settings.R
#
# T0 of the settings layout (07_documentation/agent_workflows/scope_inputs/settings_build_spec.md,
# section 8): static and unit tests of the settings, the switches and the write
# guards. Seconds; reads no data. Run from the repository root:
#   Rscript tests/check_settings.R [BASE commit, default 66c7121]
# Prints one line per check and ends with exit status 1 if any check fails.

suppressPackageStartupMessages({ library(here); library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
base <- if (length(args)) args[1] else "66c7121"
results <- list()
ok <- function(name, pass, detail = "") {
  results[[length(results) + 1]] <<- data.table(check = name, pass = isTRUE(pass), detail = detail)
  cat(sprintf("%-4s %s%s\n", if (isTRUE(pass)) "ok" else "FAIL", name, if (nzchar(detail)) paste0("  (", detail, ")") else ""))
}
env_vars <- c("WGBS_ANALYSIS", "WGBS_GENE_LEVEL_GO", "WGBS_DIAGNOSTICS")
Sys.unsetenv(env_vars)

# Settings as run_config.R gives them under the environment `env` (named values).
settings <- function(env = character(0), file = here("run_config.R")) {
  Sys.unsetenv(env_vars)
  if (length(env)) do.call(Sys.setenv, as.list(env))
  on.exit(Sys.unsetenv(env_vars))
  e <- new.env()
  source(file, local = e)
  e$config
}
stops <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")

# ---- Static ---------------------------------------------------------------------------
code <- c(list.files(c("01_pipeline", "02_R_functions", "06_diagnostics"), pattern = "\\.R$", full.names = TRUE),
          "run_config.R", "run_pipeline.R")
code <- code[!grepl("01_pipeline/optional/", code)]    # the variant step: not edited while its screen runs
hits <- unlist(lapply(code, function(f) { l <- readLines(f, warn = FALSE); f[grepl("config\\$modes|require_mode", l)] }))
ok("config$modes and require_mode() appear in no code file", length(hits) == 0, paste(unique(hits), collapse = ", "))

# ---- 16 ids x gene-level GO x diagnostics ----------------------------------------------------
source(here("02_R_functions", "analyses.R"))
tokens <- c("corrected", "cov5", "mnchisq", "ctexcl")
ids <- unique(vapply(0:15, function(k) analysis_id(tokens[bitwAnd(k, 2^(0:3)) > 0]), ""))
ok("16 analysis ids", length(ids) == 16, paste(ids, collapse = " "))
bad <- character(0)
for (id in ids) for (g in c("TRUE", "FALSE")) for (d in c("none", "all")) {
  cfg <- settings(c(WGBS_ANALYSIS = id, WGBS_GENE_LEVEL_GO = g, WGBS_DIAGNOSTICS = d))
  a <- cfg$analysis
  want_kind <- if (id == "legacy" && g == "FALSE" && d == "none") "paper" else "optional"
  fit_want <- analysis_id(setdiff(a$tokens, "ctexcl"))
  if (!(cfg$min_coverage == 10 && a$min_coverage == (if ("cov5" %in% a$tokens) 5 else 10) && a$fit_id == fit_want &&
        cfg$run_kind == want_kind && identical(cfg$addons$gene_level_go, g == "TRUE") &&
        identical(length(cfg$addons$diagnostics) > 0, d == "all"))) bad <- c(bad, paste(id, g, d))
  # paths of the analysis: never in a paper folder for an optional analysis, never under not_adopted/ for the paper in a paper run
  config <- cfg
  p <- tryCatch(analysis_input(a, "05_differential_methylation", "tables", "x.csv"), error = function(e) "error")
  if (!a$paper && !grepl("not_adopted", p)) bad <- c(bad, paste("path", id))
  if (a$paper && cfg$run_kind == "paper" && grepl("not_adopted", p)) bad <- c(bad, paste("paper path", id))
}
ok("64 combinations: min_coverage 10, analysis min_coverage, fit_id, run kind, add-ons, paths", length(bad) == 0, paste(bad, collapse = "; "))
ok("fit_id of ctexcl is legacy; of cov5-mnchisq-ctexcl is cov5-mnchisq",
   analysis_spec("ctexcl")$fit_id == "legacy" && analysis_spec("cov5-mnchisq-ctexcl")$fit_id == "cov5-mnchisq")

# ---- Invalid values stop -------------------------------------------------------------------
for (v in c("cov5-corrected", "legacy-cov5", "Cov5", "5x", " ")) {
  ok(sprintf("WGBS_ANALYSIS='%s' stops", v), stops(settings(c(WGBS_ANALYSIS = v))))
}
edited <- function(from, to) {
  f <- tempfile(fileext = ".R"); l <- readLines(here("run_config.R"))
  i <- grep(from, l, fixed = TRUE); stopifnot(length(i) == 1); l[i] <- sub(from, to, l[i], fixed = TRUE); writeLines(l, f); f
}
ok("coverage_standard = 7 stops", stops(settings(file = edited("coverage_standard          = 10L", "coverage_standard          = 7L"))))
ok("min_coverage = 5 stops (frozen)", stops(settings(file = edited("min_coverage              = 10,", "min_coverage              = 5,"))))
ok("dml_difference = 15 stops (frozen)", stops(settings(file = edited("dml_difference            = 20,", "dml_difference            = 15,"))))
ok("WGBS_GENE_LEVEL_GO=yes stops", stops(settings(c(WGBS_GENE_LEVEL_GO = "yes"))))
ok("WGBS_DIAGNOSTICS=nosuchscript stops", stops(settings(c(WGBS_DIAGNOSTICS = "nosuchscript"))))
rscript <- file.path(R.home("bin"), "Rscript")
out <- suppressWarnings(system2(rscript, c("run_pipeline.R", "--plan"), stdout = TRUE, stderr = TRUE, env = "WGBS_ANALYSIS=cov5"))
ok("run_pipeline.R refuses to start with WGBS_ANALYSIS set", !is.null(attr(out, "status")) && any(grepl("are set in this shell", out)))

# ---- Write guards ---------------------------------------------------------------------------
config <- settings(c(WGBS_ANALYSIS = "cov5"))
source(here("02_R_functions", "paths.R"))
ok("optional run: output_dir(05, tables) stops", stops(output_dir("05_differential_methylation", "tables")))
ok("optional run: manuscript_path(figures) stops", stops(manuscript_path("figures")))
ok("optional run: output_dir(05, cache) allowed", !stops(output_dir("05_differential_methylation", "cache")))
ok("optional run: manuscript_path(not_adopted, ...) allowed", !stops(manuscript_path("not_adopted", "x.csv")))
config <- settings()
ok("paper run: analysis_dir of another analysis stops", stops(analysis_dir(analysis_spec("cov5"), "04_cpg_coverage", "tables")))
rm(config)

# ---- --plan of the examples ------------------------------------------------------------------
plans <- list(
  list(args = character(0), id = "legacy", steps = "01 02 03 04 05 06 07 08 09 11 12 13"),
  list(args = "--coverage-filter", id = "corrected", steps = "04 05 06 07 08 09 13"),
  list(args = c("--coverage", "5"), id = "cov5", steps = "04 05 06 07 08 09 13"),
  list(args = "--overdispersion", id = "mnchisq", steps = "04 05 06 07 08 09 13"),
  list(args = "--exclude-ct-variants", id = "ctexcl", steps = "04 05 06 07 08 09 13"),
  list(args = "--gene-level-go", id = "legacy", steps = "10 11 13"),
  list(args = "--diagnostics", id = "legacy", steps = "13"),
  list(args = c("--coverage-filter", "--coverage", "5", "--exclude-ct-variants"), id = "corrected-cov5-ctexcl", steps = "04 05 06 07 08 09 13"),
  list(args = c("--coverage", "5", "--overdispersion", "--exclude-ct-variants"), id = "cov5-mnchisq-ctexcl", steps = "04 05 06 07 08 09 13"),
  list(args = c("--coverage", "5", "--overdispersion"), id = "cov5-mnchisq", steps = "04 05 06 07 08 09 13"),
  list(args = c("--overdispersion", "--gene-level-go"), id = "mnchisq", steps = "04 05 06 07 08 09 10 11 13"))
for (p in plans) {
  out <- suppressWarnings(system2(rscript, c("run_pipeline.R", "--plan", p$args), stdout = TRUE, stderr = TRUE))
  got_id <- sub("^  Analysis: ([^ ]+).*", "\\1", grep("^  Analysis:", out, value = TRUE))
  got_steps <- sub("^  Steps to run: ", "", grep("^  Steps to run:", out, value = TRUE))
  ok(sprintf("--plan %s", paste(p$args, collapse = " ")), identical(got_id, p$id) && identical(got_steps, p$steps),
     paste(got_id, "|", got_steps))
}
out <- suppressWarnings(system2(rscript, c("run_pipeline.R", "--coverage", "7"), stdout = TRUE, stderr = TRUE))
ok("--coverage 7 stops before a step starts", !is.null(attr(out, "status")))
out <- suppressWarnings(system2(rscript, c("run_pipeline.R", "--overdispersion", "--overdispersion", "--plan"), stdout = TRUE, stderr = TRUE))
ok("a switch given twice stops", !is.null(attr(out, "status")))

# ---- Maps -------------------------------------------------------------------------------------
read_csv_file <- function(f) utils::read.csv(f, stringsAsFactors = FALSE, colClasses = "character", na.strings = NULL, check.names = FALSE,
                                              encoding = "UTF-8")
base_file <- tempfile(fileext = ".csv")
system2("git", c("show", paste0(base, ":03_input_files/manuscript/manuscript_values_map.csv")), stdout = base_file)
base_map <- read_csv_file(base_file)
map <- read_csv_file("03_input_files/manuscript/manuscript_values_map.csv")
same <- vapply(setdiff(names(map), "expression_optional"), function(n) identical(map[[n]], base_map[[n]]), logical(1))
ok("values map: every column except expression_optional equals BASE", all(same) && identical(setdiff(names(base_map), names(map)), "expression_corrected"))
e <- map$expression_optional
e_wo_diag <- gsub('diagnostic\\(\\s*"[^"]*"\\)', "", e)
ok("values map: 263 expression_optional; no _corrected table name, no literal \"corrected\"",
   sum(nzchar(e)) == 263 && !any(grepl("_corrected", e_wo_diag)) && !any(grepl('"corrected"', e, fixed = TRUE)))
for (f in c("manuscript_values_map.csv", "manuscript_items_map.csv")) {
  raw <- readBin(file.path("03_input_files/manuscript", f), "raw", 1e7)
  n_lf <- sum(raw == as.raw(10)); n_crlf <- sum(raw[-1] == as.raw(10) & raw[-length(raw)] == as.raw(13))
  ok(sprintf("%s keeps CRLF on every line", f), n_lf > 0 && n_lf == n_crlf)
}
items <- read_csv_file("03_input_files/manuscript/manuscript_items_map.csv")
ok("items map: analyses legacy, optional, comparisons or empty",
   all(items$analysis %in% c("", "legacy", "optional", "comparisons")) && !anyDuplicated(items$destination))

res <- rbindlist(results)
cat(sprintf("\nT0: %d of %d checks pass\n", sum(res$pass), nrow(res)))
quit(status = if (all(res$pass)) 0 else 1)
