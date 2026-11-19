# 13_manuscript_values.R
#
# Purpose   Tie every number in the manuscript draft to the pipeline output that
#           regenerates it, and assemble 05_manuscript/.
#             1. Evaluate 03_input_files/manuscript/manuscript_values_map.csv.
#                Each row is one quantitative statement of the draft, with an R
#                expression that computes the same quantity from the pipeline
#                outputs as the draft's analysis did (expression_legacy) and,
#                where the audit recommends a different analysis, an expression
#                for that one (expression_corrected).
#             2. Compare the regenerated values with the manuscript and write the
#                result, one row per statement, to 05_manuscript/manuscript_values.csv.
#             3. Copy the figures and tables of the pipeline into 05_manuscript/
#                under the numbers the manuscript uses, as listed in
#                03_input_files/manuscript/manuscript_items_map.csv.
#
# Inputs    03_input_files/manuscript/manuscript_values_map.csv
#           03_input_files/manuscript/manuscript_items_map.csv
#             (item, analysis, source, destination, rows)
#           04_output/<step>/tables and figures of steps 02 to 12
# Outputs   05_manuscript/manuscript_values.csv
#           05_manuscript/manuscript_values_summary.csv
#           05_manuscript/<destination> for every row of the items map
#           04_output/13_manuscript_values/tables/expression_errors.csv (empty when all expressions run)
#
# Manuscript items: all.
# Notes     - Status of a statement:
#               matches          the regenerated value equals the manuscript value within the tolerance of the row
#               differs          it does not
#               not regenerated  the row has no expression; its note says why (for example a value
#                                that depends on information only the authors have)
#           - A row can match and still be wrong science: "matches" means the
#             pipeline reproduces the number, not that the number means what the
#             text says. The issue column points to the entry of issues.md that
#             covers it, and corrected_differs marks rows whose value changes
#             under the recommended analysis.
#           - Expressions are R code evaluated with these helpers in scope:
#               tab(step, file)    a table from 04_output/<step>/tables/ as a data.table (cached)
#               diagnostic(file)   a table from 06_diagnostics/
#               config             the settings list of run_config.R
#             An expression must return one value. A number is compared with the
#             manuscript value within the tolerance of the row; a character
#             result is compared as text; TRUE or FALSE is the answer itself, for
#             statements that are an inequality or a ranking (p < 0.05; "the
#             largest share was introns"), and TRUE counts as a match.
#           - The step stops with an error after writing its outputs if any
#             expression fails to run, so a broken row cannot pass unnoticed.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table"))

step <- "13_manuscript_values"
t0 <- start_step(step, "manuscript values and manuscript folder")

# =============================================================================
# 1. Evaluate the map
# =============================================================================
# read.csv, not fread: the expressions contain quoted strings, and base R undoes
# the doubling of quotes inside a quoted CSV field, which fread leaves in place.
map <- as.data.table(utils::read.csv(
  input_path("manuscript", "manuscript_values_map.csv"), stringsAsFactors = FALSE, na.strings = "", encoding = "UTF-8",
  colClasses = c(id = "character", section = "character", quote = "character", item = "character",
                 manuscript_value = "character", manuscript_numeric = "numeric", kind = "character", tolerance = "numeric",
                 expression_legacy = "character", expression_corrected = "character", issue = "character", note = "character")))
stopifnot(!anyDuplicated(map$id), !anyNA(map$id), !anyNA(map$section), !anyNA(map$manuscript_value))

table_cache <- new.env()
tab <- function(step_name, file) {
  key <- file.path(step_name, file)
  if (is.null(table_cache[[key]])) {
    path <- here("04_output", step_name, "tables", file)
    if (!file.exists(path)) stop("No such table: ", key)
    table_cache[[key]] <- fread(path)
  }
  copy(table_cache[[key]])
}
diagnostic <- function(file) fread(here("06_diagnostics", file))
scope <- list2env(list(tab = tab, diagnostic = diagnostic, config = config), parent = globalenv())

# Evaluate one expression; returns a list(value, error).
evaluate <- function(expression_text) {
  if (is.na(expression_text)) return(list(value = NA, error = NA_character_))
  out <- tryCatch(eval(parse(text = expression_text), envir = new.env(parent = scope)), error = function(e) e)
  if (inherits(out, "error")) return(list(value = NA, error = conditionMessage(out)))
  if (length(out) != 1) return(list(value = NA, error = sprintf("expression returned %d values, expected 1", length(out))))
  list(value = out, error = NA_character_)
}
legacy <- lapply(map$expression_legacy, evaluate)
corrected <- lapply(map$expression_corrected, evaluate)

as_text <- function(v) if (is.null(v) || is.na(v)) NA_character_ else if (is.numeric(v)) format(signif(v, 10), scientific = FALSE, trim = TRUE) else as.character(v)
agrees <- function(manuscript_numeric, manuscript_value, value, tolerance) {
  if (is.na(value)) return(NA)
  if (is.logical(value)) return(value)
  if (is.character(value)) return(identical(trimws(value), trimws(manuscript_value)))
  if (is.na(manuscript_numeric)) return(NA)
  abs(as.numeric(value) - manuscript_numeric) <= tolerance + 1e-9
}

values <- copy(map)
values[, tolerance := fifelse(is.na(tolerance), 0, tolerance)]
values[, regenerated_legacy := vapply(legacy, function(x) as_text(x$value), character(1))]
values[, regenerated_corrected := vapply(corrected, function(x) as_text(x$value), character(1))]
legacy_ok <- mapply(agrees, values$manuscript_numeric, values$manuscript_value, lapply(legacy, `[[`, "value"), values$tolerance)
corrected_ok <- mapply(agrees, values$manuscript_numeric, values$manuscript_value, lapply(corrected, `[[`, "value"), values$tolerance)
values[, status := fifelse(is.na(expression_legacy), "not regenerated",
                   fifelse(is.na(legacy_ok), "not comparable", fifelse(legacy_ok, "matches", "differs")))]
values[, corrected_differs := fifelse(is.na(expression_corrected), NA, !corrected_ok)]
values[, expression_error := fcoalesce(vapply(legacy, `[[`, character(1), "error"), vapply(corrected, `[[`, character(1), "error"))]

setcolorder(values, c("id", "section", "item", "manuscript_value", "regenerated_legacy", "status",
                      "regenerated_corrected", "corrected_differs", "issue", "note", "quote", "kind",
                      "manuscript_numeric", "tolerance", "expression_legacy", "expression_corrected", "expression_error"))
dir.create(manuscript_path(), showWarnings = FALSE)
write_table(values[, !"expression_error"], manuscript_path("manuscript_values.csv"))
summary_tab <- values[, .(statements = .N, matches = sum(status == "matches"), differs = sum(status == "differs"),
                          not_regenerated = sum(status == "not regenerated"), not_comparable = sum(status == "not comparable"),
                          change_under_recommended_analysis = sum(corrected_differs, na.rm = TRUE)), by = section]
summary_tab <- rbind(summary_tab, summary_tab[, c(list(section = "all"), lapply(.SD, sum)), .SDcols = -"section"])
write_table(summary_tab, manuscript_path("manuscript_values_summary.csv"))
print(summary_tab)

errors <- values[!is.na(expression_error), .(id, section, item, expression_legacy, expression_corrected, expression_error)]
write_table(errors, output_path(step, "tables", "expression_errors.csv"))

# =============================================================================
# 2. Figures, tables and lists under the manuscript's numbering
# =============================================================================
# A row with a `rows` expression is a CSV of which only the rows that satisfy the
# expression are written (for example the terms with p below the cutoff); all
# other rows are copied as files.
items <- as.data.table(utils::read.csv(input_path("manuscript", "manuscript_items_map.csv"), stringsAsFactors = FALSE,
                                       na.strings = "", encoding = "UTF-8", colClasses = "character"))
stopifnot(!anyDuplicated(items$destination))
items[, source_exists := file.exists(here(source))]
for (i in which(items$source_exists)) {
  dest <- manuscript_path(items$destination[i])
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (is.na(items$rows[i])) {
    file.copy(here(items$source[i]), dest, overwrite = TRUE)
  } else {
    full <- fread(here(items$source[i]))
    write_table(full[eval(parse(text = items$rows[i]))], dest)
  }
}
write_table(items, output_path(step, "tables", "manuscript_items_copied.csv"))
log_msg("%d of %d manuscript items copied into 05_manuscript/", sum(items$source_exists), nrow(items))
if (any(!items$source_exists)) {
  log_msg("not yet produced by the pipeline: %s", paste(items[source_exists == FALSE, destination], collapse = ", "))
}

finish_step(step, t0)
if (nrow(errors) > 0) {
  print(errors[, .(id, expression_error)])
  stop(nrow(errors), " expression(s) in manuscript_values_map.csv failed to run; see 04_output/13_manuscript_values/tables/expression_errors.csv")
}
