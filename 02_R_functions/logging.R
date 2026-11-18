# logging.R
# Console messages and run records.

log_msg <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), " ", sprintf(...), "\n", sep = "")
  invisible(NULL)
}

# Announce a step and return its start time.
start_step <- function(step, title) {
  cat(strrep("=", 78), "\n", step, ": ", title, "\n", strrep("=", 78), "\n", sep = "")
  Sys.time()
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

# Write a data frame as CSV; gzip when the name ends in .gz.
write_table <- function(x, path) {
  data.table::fwrite(x, path)
  invisible(path)
}
