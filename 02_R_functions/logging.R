# logging.R
# Console messages and run records.

log_msg <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), " ", sprintf(...), "\n", sep = "")
  invisible(NULL)
}

# Announce a step and return its start time. The line under the banner names the
# analysis modes of run_config.R, so that every log says which modes the run was
# configured for, also the log of a step whose results do not depend on the mode.
start_step <- function(step, title) {
  cat(strrep("=", 78), "\n", step, ": ", title, "\n", strrep("=", 78), "\n", sep = "")
  log_msg("analysis modes configured in run_config.R: %s", paste(config$modes, collapse = ", "))
  Sys.time()
}

# End a step that exists for one analysis mode when that mode is not configured.
# `step_does` says in a few words what the step does for that mode; it opens the
# message. run_pipeline.R does not start such a step without its mode; this is for
# a step started alone. Under Rscript the step ends here with exit status 0,
# before it reads or writes anything, so that a job script that runs steps one by
# one goes on. A script sourced in an interactive session cannot end early
# without an error (quit() would close the session), so there it stops with one.
require_mode <- function(mode, step_does) {
  if (mode %in% config$modes) return(invisible(TRUE))
  not_run <- sprintf("%s, and \"%s\" is not in config$modes (%s): nothing is run and no file is written",
                     step_does, mode, paste(config$modes, collapse = ", "))
  if (interactive()) stop(not_run, call. = FALSE)
  log_msg("%s", not_run)
  quit(save = "no", status = 0)
}

# Write R and package versions of the current session to the step's folder and
# report the elapsed time.
finish_step <- function(step, t0) {
  f <- file.path(output_dir(step), "session_info.txt")
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
