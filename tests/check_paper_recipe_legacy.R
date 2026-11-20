# tests/check_paper_recipe_legacy.R
#
# The paper recipe of the values map (column expression_paper_recipe of
# 03_input_files/manuscript/manuscript_values_map.csv) is expression_legacy written
# as a template, with the analysis id A where the paper's expression names the
# paper's analysis. Applied to the paper's analysis (A = "legacy", the tables of
# 04_output/<step>/) every template must give exactly the value of the paper's own
# expression, as written in 05_manuscript/manuscript_values.csv. This is the check
# that a difference between an optional analysis and the paper comes from the
# option, not from the template.
#
# Run from the repository root, after the paper run of step 13 (seconds; reads,
# writes nothing):
#   Rscript tests/check_paper_recipe_legacy.R
# Prints one line per check and ends with exit status 1 if any check fails.
# With python3 on the path it also checks that the stored column is what
# 07_documentation/toggles_migration/paper_recipe_column.py gives.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table"))
results <- list()
ok <- function(name, pass) { results[[length(results) + 1]] <<- isTRUE(pass); cat(sprintf("%-4s %s\n", if (isTRUE(pass)) "ok" else "FAIL", name)) }

map <- as.data.table(utils::read.csv(input_path("manuscript", "manuscript_values_map.csv"), stringsAsFactors = FALSE, na.strings = "",
                                     encoding = "UTF-8", colClasses = "character"))
paper <- as.data.table(utils::read.csv(manuscript_path("manuscript_values.csv"), stringsAsFactors = FALSE, na.strings = "",
                                       encoding = "UTF-8", colClasses = "character"))
ok("the map has the column expression_paper_recipe", "expression_paper_recipe" %in% names(map))
ok("a paper recipe exactly where the paper has an expression",
   identical(is.na(map$expression_paper_recipe), is.na(map$expression_legacy)))
ok("manuscript_values.csv lists the rows of the map in its order", identical(paper$id, map$id))

# The scope of the paper part of step 13, with A = "legacy".
cache <- new.env()
tab <- function(step_name, file) {
  key <- file.path(step_name, file)
  if (is.null(cache[[key]])) {
    path <- paper_path(step_name, "tables", file)
    if (!file.exists(path)) stop("No such table: ", key)
    cache[[key]] <- fread(path)
  }
  copy(cache[[key]])
}
diagnostic <- function(file) fread(here("06_diagnostics", file))
scope <- list2env(list(tab = tab, paper_tab = tab, diagnostic = diagnostic, config = config, A = "legacy"), parent = globalenv())
as_text <- function(v) if (is.null(v) || is.na(v)) NA_character_ else if (is.numeric(v)) format(signif(v, 10), scientific = FALSE, trim = TRUE) else as.character(v)
value <- vapply(map$expression_paper_recipe, function(e) {
  if (is.na(e)) return(NA_character_)
  v <- tryCatch(eval(parse(text = e), envir = new.env(parent = scope)), error = function(err) structure(conditionMessage(err), class = "failed"))
  if (inherits(v, "failed")) return(paste("error:", v))
  if (length(v) != 1) return(sprintf("error: %d values", length(v)))
  as_text(v)
}, character(1), USE.NAMES = FALSE)
with_expr <- !is.na(map$expression_paper_recipe)
failed <- map$id[with_expr & startsWith(fcoalesce(value, ""), "error:")]
ok(sprintf("all %d templates run with A = \"legacy\"%s", sum(with_expr),
           if (length(failed)) paste0(" (failed: ", paste(head(failed, 10), collapse = ", "), ")") else ""), length(failed) == 0)
differ <- map$id[with_expr & !mapply(identical, value, paper$regenerated_legacy)]
ok(sprintf("with A = \"legacy\" every template gives the value of manuscript_values.csv (%d of %d identical)%s",
           sum(with_expr) - length(differ), sum(with_expr),
           if (length(differ)) paste0("; differ: ", paste(head(differ, 10), collapse = ", ")) else ""), length(differ) == 0)

python <- Sys.which("python3")
if (nzchar(python)) {
  status <- system2(python, c("-I", shQuote(here("07_documentation", "toggles_migration", "paper_recipe_column.py")), "--check"),
                    stdout = FALSE, stderr = FALSE)
  ok("the stored column is what paper_recipe_column.py gives (--check)", status == 0)
}
cat(sprintf("\nPaper recipe, A = legacy: %d of %d checks pass\n", sum(unlist(results)), length(results)))
quit(status = if (all(unlist(results))) 0 else 1)
