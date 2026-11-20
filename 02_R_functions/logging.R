# logging.R
# Console messages and run records.

log_msg <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), " ", sprintf(...), "\n", sep = "")
  invisible(NULL)
}

# Announce a step and return its start time. The line under the banner names the
# run kind and the analysis of the run (run_config.R and the switches of
# run_pipeline.R) and, in an optional run, the folder of the analysis
# (04_output/not_adopted/<id>/), so that every log says what the run was
# configured for, also the log of a step whose results do not depend on the analysis.
start_step <- function(step, title) {
  cat(strrep("=", 78), "\n", step, ": ", title, "\n", strrep("=", 78), "\n", sep = "")
  log_msg("%s; analysis '%s'%s", if (identical(config$run_kind, "paper")) "paper run" else "OPTIONAL run (not adopted)",
          config$analysis$id,
          if (identical(config$run_kind, "paper")) " (the paper)" else paste0(", outputs under ",
            if (step == "13_manuscript_values") "04_output/not_adopted/ and 05_manuscript/not_adopted/"
            else if (grepl("^[0-9]{2}_", step)) paste0("04_output/not_adopted/", config$analysis$id, "/")
            else "06_diagnostics/",
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

# Peak memory of this R process and every process it has started (the forked
# workers of a model fit), for the log. memory_watch_start() forks a small R
# process that, once a second, adds up the proportional set size (Pss of
# /proc/<pid>/smaps_rollup) over this process and its descendants, itself
# included: the proportional set size divides a page shared after a fork among
# the processes that share it, so the sum is the memory the processes hold
# together. memory_watch_stop() returns the largest sum seen, in GB (2^30 bytes),
# or NA where it cannot be measured (not Linux, or no smaps_rollup). Sampling
# once a second misses a peak shorter than that.
memory_watch_start <- function(interval_s = 1) {
  if (.Platform$OS.type != "unix" || !file.exists("/proc/self/smaps_rollup")) return(NULL)
  stop_file <- tempfile("memory_watch_")
  root <- Sys.getpid()
  parent_of <- function(p) {
    s <- tryCatch(readLines(sprintf("/proc/%s/stat", p), warn = FALSE), error = function(e) character(0))
    if (length(s) == 0) return(NA_integer_)
    as.integer(strsplit(sub("^.*\\) ", "", s[1]), " ", fixed = TRUE)[[1]][2])   # the field after the state
  }
  tree_gb <- function() {
    pids <- as.integer(dir("/proc", pattern = "^[0-9]+$"))
    ppid <- vapply(pids, parent_of, integer(1))
    tree <- root
    repeat { add <- setdiff(pids[ppid %in% tree], tree); if (length(add) == 0) break; tree <- c(tree, add) }
    kb <- vapply(tree, function(p) {
      l <- tryCatch(readLines(sprintf("/proc/%d/smaps_rollup", p), warn = FALSE), error = function(e) character(0))
      l <- l[startsWith(l, "Pss:")]
      if (length(l) > 0) as.numeric(gsub("[^0-9]", "", l[1])) else 0
    }, numeric(1))
    sum(kb) / 1024^2
  }
  job <- parallel::mcparallel({
    peak <- 0
    repeat {
      peak <- max(peak, tryCatch(tree_gb(), error = function(e) 0))
      if (file.exists(stop_file) || !identical(parent_of("self"), root)) break   # stopped, or the R process is gone
      Sys.sleep(interval_s)
    }
    peak
  }, silent = TRUE)
  list(job = job, stop_file = stop_file)
}
memory_watch_stop <- function(watch) {
  if (is.null(watch)) return(NA_real_)
  file.create(watch$stop_file)
  on.exit(unlink(watch$stop_file))
  peak <- tryCatch(parallel::mccollect(watch$job, wait = TRUE)[[1]], error = function(e) NA_real_)
  if (is.numeric(peak) && length(peak) == 1 && peak > 0) peak else NA_real_
}
