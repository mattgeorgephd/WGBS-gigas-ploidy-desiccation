# Evaluate candidate expressions for manuscript_values_map.csv exactly as step 13 does.
# Usage: Rscript eval_expressions.R <file with one R expression per line>
# Prints, per line: the expression, its value (or the error), and its class.
suppressPackageStartupMessages({library(here); library(data.table)})
setwd("/home/user/WGBS-gigas-ploidy-desiccation")
source(here("run_config.R"))
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
lines <- readLines(commandArgs(trailingOnly = TRUE)[1])
for (i in seq_along(lines)) {
  if (!nzchar(trimws(lines[i]))) next
  out <- tryCatch(eval(parse(text = lines[i]), envir = new.env(parent = scope)), error = function(e) e)
  if (inherits(out, "error")) cat(sprintf("[%d] ERROR: %s\n     %s\n", i, conditionMessage(out), lines[i]))
  else if (length(out) != 1) cat(sprintf("[%d] ERROR: returned %d values, need exactly 1\n     %s\n", i, length(out), lines[i]))
  else cat(sprintf("[%d] %s  (%s)\n     %s\n", i, format(out, digits = 10), class(out)[1], lines[i]))
}
