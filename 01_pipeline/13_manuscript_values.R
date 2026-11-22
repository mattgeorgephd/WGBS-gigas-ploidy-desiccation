# 13_manuscript_values.R
#
# Purpose   Tie every number in the manuscript draft to the pipeline output that
#           regenerates it, and assemble 05_manuscript/.
#             1. Evaluate 03_input_files/manuscript/manuscript_values_map.csv.
#                Each row is one quantitative statement of the draft, with an R
#                expression that computes the same quantity from the pipeline
#                outputs as the draft's analysis did (expression_legacy), the same
#                computation written as a template for any analysis
#                (expression_paper_recipe; the analysis id is A) and, where the
#                audit computes the quantity another way, a template for that
#                (expression_optional, the corrected recipe).
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
#             05_manuscript/not_adopted/manuscript_values_summary_by_analysis.csv (one row per analysis and recipe)
#             05_manuscript/not_adopted/<id>/... and not_adopted/comparisons/... (items)
#             04_output/not_adopted/13_manuscript_values/tables/expression_errors_by_analysis.csv (analysis, recipe, id,
#               state "error" or "not evaluable: too few DML", expression, expression_error: the message),
#               manuscript_items_copied.csv, and session_info.txt
#
# Manuscript items: all.
# Notes     - Run kinds. The paper run evaluates expression_legacy on the paper's
#             outputs and copies the items of the paper; it never reads
#             04_output/not_adopted/. An optional run evaluates two recipes for
#             every analysis folder under 04_output/not_adopted/ that holds a
#             stamp (analysis.csv) and for the eight analyses of a single option
#             (corrected, cov5, mnchisq, ctrl, inter, ctexcl, d50, q05), whether run or not, and copies
#             their items. It does not evaluate the paper's analysis ("legacy"):
#             its values are those of 05_manuscript/manuscript_values.csv.
#           - The two recipes of an optional value (decision of the author,
#             2026-11-20, plan section 7, "Decided"):
#               paper recipe      expression_paper_recipe: the paper's own computation
#                                 (expression_legacy) on the tables of analysis A, so that a
#                                 difference from the paper's value comes only from the
#                                 option. Columns regenerated, differs_from_manuscript, status.
#                                 Generated from expression_legacy by
#                                 07_documentation/toggles_migration/paper_recipe_column.py;
#                                 with A = "legacy" every template gives the paper's value
#                                 (tests/check_paper_recipe_legacy.R). Rows that read a table
#                                 only the paper run makes (the "as drawn" clusters and the
#                                 REVIGO exports of step 11) are "not run" for every analysis.
#               corrected recipe  expression_optional: the former column expression_corrected
#                                 (the corrected mode of the audit) with A in place of
#                                 "corrected"; for some rows a corrected definition of the
#                                 quantity or another column. Columns
#                                 regenerated_corrected_recipe,
#                                 differs_from_manuscript_corrected_recipe,
#                                 status_corrected_recipe. Evaluated on the paper's outputs
#                                 these templates differ from the paper's values in 56 of 263
#                                 rows (2026-11-20). R3.1-001 and R3.1-003 take the column of
#                                 the analysis's coverage filter (cpg_in_any_library in the
#                                 legacy reading, cpg_in_any_library_unfiltered with option
#                                 corrected).
#           - Status of a statement (paper part):
#               matches          the regenerated value equals the manuscript value within the tolerance of the row
#               differs          it does not
#               not regenerated  the row has no expression; its note says why (for example a value
#                                that depends on information only the authors have)
#             In the optional part the status (of either recipe) is matches,
#             differs, not comparable (the value is NA), "not comparable: no value
#             in this analysis" (the expression gives no value: the locus, gene or
#             Table 1 row it names is not among the analysis's results), no
#             expression, "not run: analysis not run" (no folder),
#             "not run: folder stale" (its stamp describes another analysis, or
#             with ctexcl another variant list than the one on disk),
#             "not run: <step> not done (stamp)" (the stamp does not list the
#             step whose table the expression reads), "not run: <step>/<file>
#             absent" (a table the expression reads has not been made for the
#             analysis), "not evaluable: too few DML (<message>)" (the
#             expression failed and read a table of a DML set with fewer than 2
#             DML, by dml_counts_by_mode.csv of the analysis, or of a set that
#             step 09 could not test, whose Fisher table has no row, or, when
#             every set has fewer than 2 DML, any table of steps 05 to 10; or,
#             failed or not, it read the table go_gene_level_<set>_<A>.csv of a
#             set that step 10 did not test, which has no row; for example a
#             clustering of the loci in both lists at 50 points, option d50;
#             counted in column not_evaluable of the summary), or
#             "error: <message>". The optional part never stops the step: its
#             errors and its not evaluable statements, each with its message and
#             a column state, are listed in expression_errors_by_analysis.csv; a
#             warning line counts the errors and a second line the statements
#             that are not evaluable. A failing expression that reads a table of
#             a small set is not evaluable whatever made it fail, so the message
#             should be read.
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
#                                  the paper run writes; library_coverage.csv is the paper's also for an
#                                  analysis with ctexcl, whose own would differ). A table of the analysis
#                                  folder counts only while its stamp lists the step
#               paper_tab(step, file)  a table of the paper's outputs (the same as tab() in the paper part;
#                                  the paper recipe uses it only for the GO terms named in Results 3.4)
#               diagnostic(file)   a table from 06_diagnostics/
#               config             the settings list of run_config.R (in the optional part with the
#                                  min_coverage, the dml_difference and the dml_qvalue of the
#                                  analysis: with option d50 config$dml_difference is 50, so the
#                                  rows that name the threshold, or apply it to group means, use
#                                  50; with option q05 config$dml_qvalue is 0.05, so the two rows
#                                  that name the q-value cutoff, M2.5-005 and R3.2-007, give 0.05)
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
                 expression_legacy = "character", expression_paper_recipe = "character", expression_optional = "character",
                 issue = "character", note = "character")))
stopifnot(!anyDuplicated(map$id), !anyNA(map$id), !anyNA(map$section), !anyNA(map$manuscript_value))

diagnostic <- function(file) fread(here("06_diagnostics", file))

# Evaluate one expression in `scope`; returns a list(value, error, absent, empty).
# With empty_is_na (optional part) a result of length 0 is NA with empty = TRUE:
# the locus, gene or term the statement names is not in the analysis's tables.
evaluate <- function(expression_text, scope, empty_is_na = FALSE) {
  if (is.na(expression_text)) return(list(value = NA, error = NA_character_, absent = NA_character_))
  out <- tryCatch(eval(parse(text = expression_text), envir = new.env(parent = scope)),
                  table_absent = function(e) structure(list(message = conditionMessage(e)), class = "absent"),
                  error = function(e) e)
  if (inherits(out, "absent")) return(list(value = NA, error = NA_character_, absent = out$message))
  if (inherits(out, "error")) return(list(value = NA, error = conditionMessage(out), absent = NA_character_))
  if (empty_is_na && length(out) == 0) return(list(value = NA, error = NA_character_, absent = NA_character_, empty = TRUE))
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

  values <- copy(map)[, !c("expression_optional", "expression_paper_recipe")]
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
  # ---- Optional part: every analysis folder, and the eight single options -------------------
  # A table that does not depend on the analysis is read from the paper folder.
  shared_steps <- c("01_fetch_inputs", "02_alignment_qc", "03_global_methylation", "12_qpcr")
  shared_files <- c("library_coverage.csv", "genome_cpg_total.csv", "gene_go_map.csv.gz", "transcript_table.csv.gz", "go_map_id_repair.csv")
  # The paper's analysis is left out: its values are in 05_manuscript/manuscript_values.csv,
  # and the templates, applied to it, are not the paper's computations (see Notes).
  ids <- setdiff(analysis_order(c(analysis_folders(), "corrected", "cov5", "mnchisq", "ctrl", "inter", "ctexcl", "d50", "q05", config$analysis$id)), "legacy")
  ct_md5_now <- if (ct_variant_list_present()) ct_variant_list_md5() else NA_character_
  rows <- list(); error_rows <- list()
  for (A in ids) {
    a <- analysis_spec(A)
    s <- read_stamp(A)
    done <- if (is.null(s)) character(0) else strsplit(s$steps_done, ";", fixed = TRUE)[[1]]
    folder_state <- if (is.null(s)) "not run: analysis not run" else if (
      s$filter != a$filter || s$min_coverage != format(a$min_coverage) || s$overdispersion != a$overdispersion ||
      s$ct_excluded != as.character(a$ct_excluded) || s$dml_difference != format(a$dml_difference) ||
      s$dml_qvalue != format(a$dml_qvalue) || s$contrast_variant != a$contrast_variant ||
      (a$ct_excluded && !is.na(ct_md5_now) && !identical(s$ct_variant_list_md5, ct_md5_now))) "not run: folder stale" else ""
    recipes <- c(paper = "expression_paper_recipe", corrected = "expression_optional")
    if (nzchar(folder_state)) {
      res <- lapply(recipes, function(r) rep(list(list(value = NA, error = NA_character_, absent = NA_character_)), nrow(map)))
    } else {
      table_cache <- new.env()
      read_log <- new.env(); read_log$files <- character(0)      # tables of the analysis read by the expression being evaluated
      tab <- function(step_name, file) {
        key <- file.path(step_name, file)
        read_log$files <- c(read_log$files, key)
        if (is.null(table_cache[[key]])) {
          absent <- function(message) stop(structure(class = c("table_absent", "error", "condition"), list(message = message, call = NULL)))
          path <- not_adopted_path(A, step_name, "tables", file)
          if (!file.exists(path) && (step_name %in% shared_steps || file %in% shared_files)) {
            path <- paper_path(step_name, "tables", file)
          } else if (!substr(step_name, 1, 2) %in% done) {
            absent(paste0("not run: ", step_name, " not done (stamp)"))
          }
          if (!file.exists(path)) absent(paste0("not run: ", step_name, "/", file, " absent"))
          table_cache[[key]] <- fread(path)
        }
        copy(table_cache[[key]])
      }
      paper_tab <- function(step_name, file) {
        path <- paper_path(step_name, "tables", file)
        if (!file.exists(path)) stop("No such table of the paper: ", file.path(step_name, file))
        fread(path)
      }
      # DML sets too small for some statements: fewer than 2 DML (dml_counts_by_mode.csv of the
      # analysis; 2 is the fewest that can be clustered or compared), or, for a table of step 09,
      # a set that step 09 could not test (its Fisher table has no row: no DML in a gene of the
      # legacy universe with a term of the ontology). An expression that fails and read a table
      # of such a set (named by the set: dml_ploidy_<A>.csv, legacy_fisher_shared_BP.csv, ...),
      # or, when every set has fewer than 2 DML, any table of steps 05 to 10 of the analysis, is
      # "not evaluable: too few DML (<message>)", not an error; the message is kept.
      counts_file <- not_adopted_path(A, "05_differential_methylation", "tables", "dml_counts_by_mode.csv")
      set_n <- if (file.exists(counts_file) && "05" %in% done) {
        k <- fread(counts_file)[mode == A]
        if (nrow(k) == 1) c(ploidy = k$ploidy_dml, heat = k$heat_dml, shared = k$shared) else NULL
      }
      reads_small_set <- function(keys) {
        if (is.null(set_n) || length(keys) == 0) return(FALSE)
        file <- basename(keys)
        # the DML set of a table: dml_ploidy_<A>.csv, legacy_fisher_shared_BP.csv, genes_with_dml_heat_<A>.csv, ...
        set_of <- vapply(regmatches(file, regexec("_(ploidy|heat|shared)[_.]", file)), function(m) if (length(m)) m[2] else NA_character_, "")
        of_set <- !is.na(set_of)
        small <- of_set & unname(set_n[set_of]) %in% c(0, 1)
        all_small <- !of_set & all(set_n < 2) & grepl("^(0[5-9]|10)_", keys)
        untested <- of_set & startsWith(file, "legacy_fisher_") & vapply(keys, function(k) {
          x <- table_cache[[k]]; !is.null(x) && nrow(x) == 0 }, logical(1))
        any(small | untested | all_small)
      }
      # Sets that step 10 did not test (no DML gene carries a term of the ontology): their table
      # go_gene_level_<set>_<A>.csv has no row. Any value computed from such a table, failed or
      # not, says nothing about the set, and is "not evaluable" as well.
      sets_not_tested_by_step10 <- function(keys) {
        k <- keys[grepl("^10_go_enrichment_gene_level/go_gene_level_(ploidy|heat|shared)_", keys)]
        empty <- vapply(k, function(x) { t <- table_cache[[x]]; !is.null(t) && nrow(t) == 0 }, logical(1))
        unique(sub("^10_go_enrichment_gene_level/go_gene_level_(ploidy|heat|shared)_.*", "\\1", k[empty]))
      }
      cfg_a <- modifyList(config, list(min_coverage = a$min_coverage, dml_difference = a$dml_difference, dml_qvalue = a$dml_qvalue))
      scope <- list2env(list(tab = tab, paper_tab = paper_tab, diagnostic = diagnostic, config = cfg_a, A = A), parent = globalenv())
      res <- lapply(recipes, function(r) lapply(map[[r]], function(e) {
        read_log$files <- character(0)
        x <- evaluate(e, scope, empty_is_na = TRUE)
        # The reason a statement is not evaluable, or NA. The message of a failed expression is kept.
        x$not_evaluable <- NA_character_
        if (!is.na(x$error) && reads_small_set(read_log$files)) {
          x$not_evaluable <- x$error
        } else if (is.na(x$absent)) {
          untested <- sets_not_tested_by_step10(read_log$files)
          if (length(untested) > 0) {
            x$not_evaluable <- sprintf("step 10 did not test the %s set: no DML gene carries a term of the ontology%s",
                                       paste(untested, collapse = " and "), if (is.na(x$error)) "" else paste0("; ", x$error))
          }
        }
        x
      }))
    }
    # Value, whether it differs from the manuscript, and status of each row under one recipe.
    judge <- function(r, expressions) {
      value <- lapply(r, `[[`, "value")
      field <- function(name) vapply(r, function(x) if (is.null(x[[name]])) NA_character_ else x[[name]], character(1))
      empty <- vapply(r, function(x) isTRUE(x$empty), logical(1))
      not_evaluable <- vapply(r, function(x) if (is.null(x$not_evaluable)) NA_character_ else x$not_evaluable, character(1))
      too_few <- !is.na(not_evaluable)
      ok <- mapply(agrees, map$manuscript_numeric, map$manuscript_value, value, map$tolerance)
      status <- if (nzchar(folder_state)) rep(folder_state, nrow(map)) else
        fifelse(is.na(expressions), "no expression",
        fifelse(!is.na(field("absent")), field("absent"),
        fifelse(too_few, paste0("not evaluable: too few DML (", not_evaluable, ")"),
        fifelse(!is.na(field("error")), paste("error:", field("error")),
        fifelse(empty, "not comparable: no value in this analysis",
        fifelse(is.na(ok), "not comparable", fifelse(ok, "matches", "differs")))))))
      list(regenerated = vapply(value, as_text, character(1)),
           differs = fifelse(status %in% c("matches", "differs"), status == "differs", NA), status = status,
           error = fifelse(too_few, not_evaluable, field("error")),
           state = fifelse(too_few, "not evaluable: too few DML", fifelse(is.na(field("error")), NA_character_, "error")))
    }
    p <- judge(res$paper, map$expression_paper_recipe)
    k <- judge(res$corrected, map$expression_optional)
    rows[[A]] <- data.table(id = map$id, section = map$section, item = map$item, manuscript_value = map$manuscript_value, analysis = A,
                            regenerated = p$regenerated, differs_from_manuscript = p$differs, status = p$status,
                            regenerated_corrected_recipe = k$regenerated, differs_from_manuscript_corrected_recipe = k$differs,
                            status_corrected_recipe = k$status)
    # Every failed expression and every statement that is not evaluable, with its message.
    err <- rbind(data.table(analysis = A, recipe = "paper", id = map$id, state = p$state, expression = map$expression_paper_recipe,
                            expression_error = p$error),
                 data.table(analysis = A, recipe = "corrected", id = map$id, state = k$state, expression = map$expression_optional,
                            expression_error = k$error))
    if (any(!is.na(err$state))) error_rows[[A]] <- err[!is.na(state)]
  }
  by_analysis <- rbindlist(rows)
  dir.create(manuscript_path("not_adopted"), recursive = TRUE, showWarnings = FALSE)
  write_table(by_analysis, manuscript_path("not_adopted", "manuscript_values_by_analysis.csv"))
  # One row per analysis and recipe, over the statements that have an expression of that recipe.
  count_recipe <- function(recipe_name, status_column, expressions) {
    x <- data.table(analysis = by_analysis$analysis, s = by_analysis[[status_column]], e = expressions[match(by_analysis$id, map$id)])
    x[!is.na(e), .(recipe = recipe_name, statements_with_expression = .N, matches = sum(s == "matches"), differs = sum(s == "differs"),
        not_comparable = sum(startsWith(s, "not comparable")), not_run = sum(startsWith(s, "not run")),
        not_evaluable = sum(startsWith(s, "not evaluable")), errors = sum(startsWith(s, "error")),
        state = if (all(startsWith(s, "not run: analysis"))) "analysis not run" else
          if (all(s == "not run: folder stale")) "folder stale" else "evaluated"), by = analysis]
  }
  summary_tab <- rbind(count_recipe("paper", "status", map$expression_paper_recipe),
                       count_recipe("corrected", "status_corrected_recipe", map$expression_optional))
  summary_tab <- summary_tab[order(match(analysis, unique(by_analysis$analysis)), match(recipe, c("paper", "corrected")))]
  write_table(summary_tab, manuscript_path("not_adopted", "manuscript_values_summary_by_analysis.csv"))
  print(summary_tab)
  errors <- if (length(error_rows)) rbindlist(error_rows) else
    data.table(analysis = character(0), recipe = character(0), id = character(0), state = character(0), expression = character(0),
               expression_error = character(0))
  write_table(errors, output_path(file.path("not_adopted", step), "tables", "expression_errors_by_analysis.csv"))
  real <- errors[state == "error"]
  if (nrow(real) > 0) {
    log_msg("WARNING: %d optional expression(s) failed to run (%s); see 04_output/not_adopted/13_manuscript_values/tables/expression_errors_by_analysis.csv",
            nrow(real), paste(sprintf("%s %s", real$analysis, real$id), collapse = ", "))
  }
  if (nrow(errors) > nrow(real)) {
    log_msg("%d optional statement(s) not evaluable because a DML set has too few DML (%s); their messages are in 04_output/not_adopted/13_manuscript_values/tables/expression_errors_by_analysis.csv (state 'not evaluable: too few DML')",
            nrow(errors) - nrow(real), paste(unique(errors[state != "error", analysis]), collapse = ", "))
  }
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
    # The folder not_adopted/legacy/ holds only what an optional run makes for the
    # paper's analysis (steps 10 and 11, gene level); its other items are the paper's.
    if (A == "legacy") x <- x[file.exists(here(source))]
    x
  })), items[analysis == "comparisons"])
}
# A source in an analysis folder counts only while the stamp of the folder lists its step
# (the comparison folder has no stamp; rebuild_comparison() reads only folders that list the step).
in_folder <- !is.na(items$analysis) & items$analysis != "comparisons" & grepl("^04_output/not_adopted/[^/]+/[0-9]{2}_", items$source)
items[, source_exists := file.exists(here(source))]
for (i in which(in_folder)) {
  set(items, i, "source_exists", items$source_exists[i] && stamp_lists_step(sub("^04_output/not_adopted/([^/]+)/.*", "\\1", items$source[i]),
                                                                              sub("^04_output/not_adopted/[^/]+/([0-9]{2})_.*", "\\1", items$source[i])))
}
# No destination may leave 05_manuscript/ (or, in an optional run, 05_manuscript/not_adopted/).
dotdot <- vapply(strsplit(items$destination, "[/\\\\]"), function(p) any(p == ".."), logical(1))
if (any(dotdot)) stop("manuscript_items_map.csv: a destination may not contain '..': ", paste(items$destination[dotdot], collapse = ", "))
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
