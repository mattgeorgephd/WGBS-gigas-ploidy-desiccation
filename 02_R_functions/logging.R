# logging.R
# Console messages and run records.

log_msg <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), " ", sprintf(...), "\n", sep = "")
  invisible(NULL)
}

# Announce a step and return its start time. The line under the banner names the
# run kind and the analysis of the run (run_config.R and the switches of
# run_pipeline.R), so that every log says what the run was configured for, also
# the log of a step whose results do not depend on the analysis.
start_step <- function(step, title) {
  cat(strrep("=", 78), "\n", step, ": ", title, "\n", strrep("=", 78), "\n", sep = "")
  log_msg("%s; analysis '%s'%s", if (identical(config$run_kind, "paper")) "paper run" else "OPTIONAL run (not adopted)",
          config$analysis$id,
          if (identical(config$run_kind, "paper")) " (the paper)" else paste0(", outputs under 04_output/not_adopted/",
            if (isTRUE(config$addons$gene_level_go)) "; gene-level GO on" else "",
            if (length(config$addons$diagnostics)) paste0("; diagnostics: ", paste(config$addons$diagnostics, collapse = ",")) else ""))
  Sys.time()
}

# Write R and package versions of the current session to the step's folder and
# report the elapsed time.
finish_step <- function(step, t0) {
  # The paper run writes to the paper folder; an optional run to the folder of its
  # analysis (04_output/not_adopted/<id>/<step>/) or, for a step that is not an
  # analysis step, to 04_output/not_adopted/<step>/.
  folder <- if (!identical(config$run_kind, "optional") || startsWith(step, "not_adopted")) output_dir(step)
            else if (step %in% c("13_manuscript_values")) output_dir(file.path("not_adopted", step))
            else analysis_dir(config$analysis, step)
  f <- file.path(folder, "session_info.txt")
  writeLines(c(paste("Step:", step), paste("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
               capture.output(print(utils::sessionInfo()))), f)
  log_msg("%s finished in %.1f min", step, as.numeric(difftime(Sys.time(), t0, units = "mins")))
  invisible(f)
}

# Write a data frame as CSV; gzip when the name ends in .gz. The bytes must not
# depend on the platform or on the options of the R session, because the tables
# are compared byte for byte after a rerun. Left to its defaults fwrite() writes
# CRLF line endings on Windows, takes the notation of numbers from the option
# scipen of the session (a personal R profile often sets it), and writes TRUE
# and FALSE as 1 and 0 under the option datatable.logical01. A gzipped table
# still differs between operating systems in one byte of the gzip header, which
# names the system; its decompressed content does not.
write_table <- function(x, path) {
  data.table::fwrite(x, path, eol = "\n", scipen = 0, logical01 = FALSE)
  invisible(path)
}
