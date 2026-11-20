# 13_manuscript_values.R
#
# Purpose   Tie every number in the manuscript draft to the pipeline output that
#           regenerates it, and assemble 05_manuscript/.
#             1. Evaluate 03_input_files/manuscript/manuscript_values_map.csv.
#                Each row is one quantitative statement of the draft, with an R
#                expression that computes the same quantity from the pipeline
#                outputs as the draft's analysis did (expression_legacy) and,
#                where an optional analysis gives the quantity another way, a
#                template for it (expression_optional; the analysis id is A).
#             2. Compare the regenerated values with the manuscript and write the
#                result, one row per statement, to 05_manuscript/manuscript_values.csv.
#             3. Copy the figures and tables of the pipeline into 05_manuscript/
#                under the numbers the manuscript uses, as listed in
#                03_input_files/manuscript/manuscript_items_map.csv.
#
# Run       Rscript 01_pipeline/13_manuscript_values.R
#           (after every other step; a few seconds. Rerun it after any step that
#           changes an output.)
#
# Inputs    03_input_files/manuscript/manuscript_values_map.csv
#           03_input_files/manuscript/manuscript_items_map.csv
#             (item, analysis, source, destination, rows)
#           04_output/<step>/tables and figures of steps 02 to 12
#           04_output/not_adopted/<id>/ (optional part)
#           the tables of 06_diagnostics/ that an expression names
# Outputs   Paper run (paper part):
#             05_manuscript/manuscript_values.csv
#             05_manuscript/manuscript_values_summary.csv
#             05_manuscript/<destination> for every row of the items map of the paper
#             04_output/13_manuscript_values/tables/expression_errors.csv (empty when all expressions run)
#             04_output/13_manuscript_values/tables/manuscript_items_copied.csv
#             04_output/13_manuscript_values/session_info.txt
#           Optional run (optional part; NOT ADOPTED analyses):
#             05_manuscript/not_adopted/manuscript_values_by_analysis.csv
#             05_manuscript/not_adopted/manuscript_values_summary_by_analysis.csv
#             05_manuscript/not_adopted/<id>/... and not_adopted/comparisons/... (items)
#             04_output/not_adopted/13_manuscript_values/tables/expression_errors_by_analysis.csv,
#               manuscript_items_copied.csv, and session_info.txt
#
# Manuscript items: all.
# Notes     - Run kinds. The paper run evaluates expression_legacy on the paper's
#             outputs and copies the items of the paper; it never reads
#             04_output/not_adopted/. An optional run evaluates
#             expression_optional for every analysis folder under
#             04_output/not_adopted/ that holds a stamp (analysis.csv) and for
#             the four analyses of a single option (corrected, cov5, mnchisq,
#             ctexcl), whether run or not, and copies their items.
#           - Status of a statement (paper part):
#               matches          the regenerated value equals the manuscript value within the tolerance of the row
#               differs          it does not
#               not regenerated  the row has no expression; its note says why (for example a value
#                                that depends on information only the authors have)
#             In the optional part the status is matches, differs, not
#             comparable, no expression, "not run: analysis not run" (no folder),
#             "not run: folder stale" (its stamp describes another analysis),
#             "not run: <step>/<file> absent" (a table the expression reads has
#             not been made for the analysis), or "error: <message>". The
#             optional part never stops the step: its errors are listed in
#             expression_errors_by_analysis.csv.
#           - A row can match and still be wrong science: "matches" means the
#             pipeline reproduces the number, not that the number means what the
#             text says. The issue column points to the entry of issues.md that
#             covers it. The optional analyses are not analyses that issues.md
#             recommends for the DML themselves (I-06) or for their location (I-41).
#           - Expressions are R code evaluated with these helpers in scope:
#               tab(step, file)    a table from 04_output/<step>/tables/ as a data.table (cached); in the
#                                  optional part from 04_output/not_adopted/<A>/<step>/tables/, or from
#                                  the paper folder for a table that does not depend on the analysis
#                                  (steps 01 to 03 and 12, and the tables of steps 04 and 08 that only
#                                  the paper run writes) and, for A = legacy, for every table
#               paper_tab(step, file)  a table of the paper's outputs (the same as tab() in the paper part)
#               diagnostic(file)   a table from 06_diagnostics/
#               config             the settings list of run_config.R (in the optional part with the
#                                  min_coverage of the analysis)
#               A                  (optional part) the analysis id
#             An expression must return one value. A number is compared with the
#             manuscript value within the tolerance of the row; a character
#             result is compared as text; TRUE or FALSE is the answer itself, for
#             statements that are an inequality or a ranking (p < 0.05; "the
#             largest share was introns"), and TRUE counts as a match.
#           - In the paper part the step stops with an error after writing its
#             outputs if any expression fails to run or reads a table that is
#             absent, so a broken row cannot pass unnoticed.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table"))

step <- "13_manuscript_values"
t0 <- start_step(step, "manuscript values and manuscript folder")
paper_run <- identical(config$run_kind, "paper")

# =============================================================================
# 1. Evaluate the map
# =============================================================================
# read.csv, not fread: the expressions contain quoted strings, and base R undoes
# the doubling of quotes inside a quoted CSV field, which fread leaves in place.
map <- as.data.table(utils::read.csv(
  input_path("manuscript", "manuscript_values_map.csv"), stringsAsFactors = FALSE, na.strings = "", encoding = "UTF-8",
  colClasses = c(id = "character", section = "character", quote = "character", item = "character",
                 manuscript_value = "character", manuscript_numeric = "numeric", kind = "character", tolerance = "numeric",
                 expression_legacy = "character", expression_optional = "character", issue = "character", note = "character")))
stopifnot(!anyDuplicated(map$id), !anyNA(map$id), !anyNA(map$section), !anyNA(map$manuscript_value))

diagnostic <- function(file) fread(here("06_diagnostics", file))

# Evaluate one expression in `scope`; returns a list(value, error, absent).
evaluate <- function(expression_text, scope) {
  if (is.na(expression_text)) return(list(value = NA, error = NA_character_, absent = NA_character_))
  out <- tryCatch(eval(parse(text = expression_text), envir = new.env(parent = scope)),
                  table_absent = function(e) structure(list(message = conditionMessage(e)), class = "absent"),
                  error = function(e) e)
  if (inherits(out, "absent")) return(list(value = NA, error = NA_character_, absent = out$message))
  if (inherits(out, "error")) return(list(value = NA, error = conditionMessage(out), absent = NA_character_))
  if (length(out) != 1) return(list(value = NA, error = sprintf("expression returned %d values, expected 1", length(out)), absent = NA_character_))
  list(value = out, error = NA_character_, absent = NA_character_)
}
as_text <- function(v) if (is.null(v) || is.na(v)) NA_character_ else if (is.numeric(v)) format(signif(v, 10), scientific = FALSE, trim = TRUE) else as.character(v)
agrees <- function(manuscript_numeric, manuscript_value, value, tolerance) {
  if (is.na(value)) return(NA)
  if (is.logical(value)) return(value)
  if (is.character(value)) return(identical(trimws(value), trimws(manuscript_value)))
  if (is.na(manuscript_numeric)) return(NA)
  abs(as.numeric(value) - manuscript_numeric) <= tolerance + 1e-9
}
map[, tolerance := fifelse(is.na(tolerance), 0, tolerance)]

if (paper_run) {
  # ---- Paper part: the paper's outputs only -------------------------------------------------
  table_cache <- new.env()
  tab <- function(step_name, file) {
    key <- file.path(step_name, file)
    if (is.null(table_cache[[key]])) {
      path <- paper_path(step_name, "tables", file)
      if (!file.exists(path)) stop("No such table: ", key)
      table_cache[[key]] <- fread(path)
    }
    copy(table_cache[[key]])
  }
  scope <- list2env(list(tab = tab, paper_tab = tab, diagnostic = diagnostic, config = config), parent = globalenv())
  legacy <- lapply(map$expression_legacy, evaluate, scope = scope)

  values <- copy(map)[, !"expression_optional"]
  values[, regenerated_legacy := vapply(legacy, function(x) as_text(x$value), character(1))]
  legacy_ok <- mapply(agrees, values$manuscript_numeric, values$manuscript_value, lapply(legacy, `[[`, "value"), values$tolerance)
  values[, status := fifelse(is.na(expression_legacy), "not regenerated",
                     fifelse(is.na(legacy_ok), "not comparable", fifelse(legacy_ok, "matches", "differs")))]
  values[, expression_error := vapply(legacy, `[[`, character(1), "error")]

  setcolorder(values, c("id", "section", "item", "manuscript_value", "regenerated_legacy", "status",
                        "issue", "note", "quote", "kind", "manuscript_numeric", "tolerance", "expression_legacy", "expression_error"))
  dir.create(manuscript_path(), showWarnings = FALSE)
  write_table(values[, !"expression_error"], manuscript_path("manuscript_values.csv"))
  summary_tab <- values[, .(statements = .N, matches = sum(status == "matches"), differs = sum(status == "differs"),
                            not_regenerated = sum(status == "not regenerated"), not_comparable = sum(status == "not comparable")),
                        by = section]
  summary_tab <- rbind(summary_tab, summary_tab[, c(list(section = "all"), lapply(.SD, sum)), .SDcols = -"section"])
  write_table(summary_tab, manuscript_path("manuscript_values_summary.csv"))
  print(summary_tab)

  errors <- values[!is.na(expression_error), .(id, section, item, expression_legacy, expression_error)]
  write_table(errors, output_path(step, "tables", "expression_errors.csv"))
} else {
  # ---- Optional part: every analysis folder, and the four single options --------------------
  # A table that does not depend on the analysis is read from the paper folder.
  shared_steps <- c("01_fetch_inputs", "02_alignment_qc", "03_global_methylation", "12_qpcr")
  shared_files <- c("library_coverage.csv", "genome_cpg_total.csv", "gene_go_map.csv.gz", "transcript_table.csv.gz", "go_map_id_repair.csv")
  ids <- analysis_order(c(analysis_folders(), "corrected", "cov5", "mnchisq", "ctexcl", config$analysis$id))
  rows <- list(); error_rows <- list()
  for (A in ids) {
    a <- analysis_spec(A)
    s <- read_stamp(A)
    folder_state <- if (is.null(s) && !(A == "legacy" && A == config$analysis$id)) "not run: analysis not run" else if (!is.null(s) && (
      s$filter != a$filter || s$min_coverage != format(a$min_coverage) || s$overdispersion != a$overdispersion ||
      s$ct_excluded != as.character(a$ct_excluded))) "not run: folder stale" else ""
    if (nzchar(folder_state)) {
      res <- rep(list(list(value = NA, error = NA_character_, absent = NA_character_)), nrow(map))
    } else {
      table_cache <- new.env()
      tab <- function(step_name, file) {
        key <- file.path(step_name, file)
        if (is.null(table_cache[[key]])) {
          path <- not_adopted_path(A, step_name, "tables", file)
          if (!file.exists(path) && (A == "legacy" || step_name %in% shared_steps || file %in% shared_files)) {
            path <- paper_path(step_name, "tables", file)
          }
          if (!file.exists(path)) stop(structure(class = c("table_absent", "error", "condition"),
                                                 list(message = paste0("not run: ", step_name, "/", file, " absent"), call = NULL)))
          table_cache[[key]] <- fread(path)
        }
        copy(table_cache[[key]])
      }
      paper_tab <- function(step_name, file) {
        path <- paper_path(step_name, "tables", file)
        if (!file.exists(path)) stop("No such table of the paper: ", file.path(step_name, file))
        fread(path)
      }
      cfg_a <- modifyList(config, list(min_coverage = a$min_coverage))
      scope <- list2env(list(tab = tab, paper_tab = paper_tab, diagnostic = diagnostic, config = cfg_a, A = A), parent = globalenv())
      res <- lapply(map$expression_optional, evaluate, scope = scope)
    }
    value <- lapply(res, `[[`, "value")
    ok <- mapply(agrees, map$manuscript_numeric, map$manuscript_value, value, map$tolerance)
    status <- if (nzchar(folder_state)) rep(folder_state, nrow(map)) else
      fifelse(is.na(map$expression_optional), "no expression",
      fifelse(!is.na(vapply(res, `[[`, character(1), "absent")), vapply(res, `[[`, character(1), "absent"),
      fifelse(!is.na(vapply(res, `[[`, character(1), "error")), paste("error:", vapply(res, `[[`, character(1), "error")),
      fifelse(is.na(ok), "not comparable", fifelse(ok, "matches", "differs")))))
    rows[[A]] <- data.table(id = map$id, section = map$section, item = map$item, manuscript_value = map$manuscript_value, analysis = A,
                            regenerated = vapply(value, as_text, character(1)),
                            differs_from_manuscript = fifelse(status %in% c("matches", "differs"), status == "differs", NA),
                            status = status)
    err <- vapply(res, `[[`, character(1), "error")
    if (any(!is.na(err))) error_rows[[A]] <- data.table(analysis = A, id = map$id, expression_optional = map$expression_optional,
                                                         expression_error = err)[!is.na(expression_error)]
  }
  by_analysis <- rbindlist(rows)
  dir.create(manuscript_path("not_adopted"), recursive = TRUE, showWarnings = FALSE)
  write_table(by_analysis, manuscript_path("not_adopted", "manuscript_values_by_analysis.csv"))
  summary_tab <- by_analysis[!is.na(map$expression_optional[match(id, map$id)]),
                             .(statements_with_expression = .N, matches = sum(status == "matches"), differs = sum(status == "differs"),
                               not_comparable = sum(status == "not comparable"), not_run = sum(startsWith(status, "not run")),
                               errors = sum(startsWith(status, "error")),
                               state = if (all(startsWith(status, "not run: analysis"))) "analysis not run" else
                                 if (all(status == "not run: folder stale")) "folder stale" else "evaluated"), by = analysis]
  write_table(summary_tab, manuscript_path("not_adopted", "manuscript_values_summary_by_analysis.csv"))
  print(summary_tab)
  errors <- if (length(error_rows)) rbindlist(error_rows) else
    data.table(analysis = character(0), id = character(0), expression_optional = character(0), expression_error = character(0))
  write_table(errors, output_path(file.path("not_adopted", step), "tables", "expression_errors_by_analysis.csv"))
}

# =============================================================================
# 2. Figures, tables and lists under the manuscript's numbering
# =============================================================================
# A row with a `rows` expression is a CSV of which only the rows that satisfy the
# expression are written (for example the terms with p below the cutoff); all
# other rows are copied as files. The paper run copies the items of the paper
# (analysis "legacy" or empty); an optional run the templates of analysis
# "optional" (one copy per analysis folder whose source exists, {A} replaced by
# its id) and the items of the comparison folder ("comparisons").
items <- as.data.table(utils::read.csv(input_path("manuscript", "manuscript_items_map.csv"), stringsAsFactors = FALSE,
                                       na.strings = "", encoding = "UTF-8", colClasses = "character"))
stopifnot(!anyDuplicated(items$destination), all(is.na(items$analysis) | items$analysis %in% c("legacy", "optional", "comparisons")))
if (paper_run) {
  items <- items[is.na(analysis) | analysis == "legacy"]
} else {
  templ <- items[analysis == "optional"]
  items <- rbind(rbindlist(lapply(analysis_folders(), function(A) {
    x <- copy(templ); x[, `:=`(analysis = A, source = gsub("{A}", A, source, fixed = TRUE), destination = gsub("{A}", A, destination, fixed = TRUE))]
  })), items[analysis == "comparisons"])
}
items[, source_exists := file.exists(here(source))]
for (i in which(items$source_exists)) {
  dest <- here("05_manuscript", items$destination[i])
  if (!paper_run) stopifnot(startsWith(items$destination[i], "not_adopted/"))
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (is.na(items$rows[i])) {
    file.copy(here(items$source[i]), dest, overwrite = TRUE)
  } else {
    full <- fread(here(items$source[i]))
    write_table(full[eval(parse(text = items$rows[i]))], dest)
  }
}
write_table(items, if (paper_run) output_path(step, "tables", "manuscript_items_copied.csv") else
                     output_path(file.path("not_adopted", step), "tables", "manuscript_items_copied.csv"))
log_msg("%d of %d manuscript items copied into 05_manuscript/%s", sum(items$source_exists), nrow(items), if (paper_run) "" else "not_adopted/")
if (any(!items$source_exists)) {
  log_msg("not produced by the pipeline (or by the analyses run): %s", paste(items[source_exists == FALSE, destination], collapse = ", "))
}

finish_step(step, t0)
if (paper_run && nrow(errors) > 0) {
  print(errors[, .(id, expression_error)])
  stop(nrow(errors), " expression(s) in manuscript_values_map.csv failed to run; see 04_output/13_manuscript_values/tables/expression_errors.csv")
}
