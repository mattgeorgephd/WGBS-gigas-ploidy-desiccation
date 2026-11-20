# run_pipeline.R
#
# Orchestrator. Runs the numbered scripts in 01_pipeline/ in order, each in a
# fresh R process, and stops at the first failure.
#
#   Rscript run_pipeline.R                    all steps
#   Rscript run_pipeline.R --steps 5:7        steps 05, 06, 07
#   Rscript run_pipeline.R --steps 3,12       steps 03 and 12
#   Rscript run_pipeline.R --steps 2:3,6:10   steps 02, 03 and 06 to 10
#   Rscript run_pipeline.R --list             show the steps and exit
#
# The selection after --steps is one word without spaces. Anything else on the
# command line stops the run before a step starts.
#
# Settings are in run_config.R. They are read and checked here before the first
# step starts, and the analysis modes of the run are printed at the start and at
# the end. The first run downloads 10.5 GB from gannet (step 01; the inputs of
# steps 04 and 05) and the 31 MB GO.db archive of step 09. Step 04 needs about
# 9.5 GB of memory and step 05 about 10 GB with WGBS_CORES=4.
#
# The repository is taken to be the folder this file lies in, so the file can be
# started from any working directory. The steps find the repository through the
# working directory (here::here()); each step is therefore started with the
# repository as working directory.

args <- commandArgs(trailingOnly = TRUE)

# ---- The repository and its steps ------------------------------------------------
# Under Rscript the path of this file is on the command line, with each space
# written as "~+~". When the file is sourced from an R session instead, the
# repository is looked for upward from the working directory, as the steps do.
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
this_file <- gsub("~+~", " ", this_file, fixed = TRUE)
root <- if (length(this_file) == 1 && file.exists(this_file)) dirname(normalizePath(this_file)) else here::here()
scripts <- sort(list.files(file.path(root, "01_pipeline"), pattern = "^[0-9]{2}_.*\\.R$", full.names = TRUE))
if (length(scripts) == 0 || !file.exists(file.path(root, "run_config.R"))) {
  stop("No pipeline found in ", root, ": it must hold run_config.R and the step scripts in 01_pipeline/. ",
       "Start run_pipeline.R with Rscript, or source it with the working directory inside the repository (it is ",
       getwd(), ").", call. = FALSE)
}
numbers <- as.integer(substr(basename(scripts), 1, 2))
if (anyDuplicated(numbers)) stop("Two scripts in 01_pipeline/ have the same step number.", call. = FALSE)

# ---- Command line ------------------------------------------------------------------
usage <- paste("Usage: Rscript run_pipeline.R [--list | --steps <selection>]\n",
               " <selection> is one word without spaces: step numbers and ranges joined by commas,",
               "for example 5:7 or 3,12 or 2:3,6:10.")
if (identical(args, "--help") || identical(args, "-h")) {
  cat(usage, "\n", sep = "")
  quit(status = 0)
}
if (identical(args, "--list")) {
  cat(sprintf("%02d  %s\n", numbers, basename(scripts)), sep = "")
  quit(status = 0)
}
selected <- numbers
if (length(args) > 0) {
  part <- "[0-9]{1,2}(:[0-9]{1,2})?"
  if (length(args) != 2 || args[1] != "--steps" || !grepl(sprintf("^%s(,%s)*$", part, part), args[2])) {
    stop("Cannot read the argument(s): ", paste(shQuote(args), collapse = " "), "\n", usage, call. = FALSE)
  }
  selected <- unlist(lapply(strsplit(args[2], ",")[[1]], function(p) {
    r <- as.integer(strsplit(p, ":")[[1]])
    if (length(r) == 2 && r[1] > r[2]) stop("A range must run upward: ", p, "\n", usage, call. = FALSE)
    if (length(r) == 2) seq(r[1], r[2]) else r
  }))
  unknown <- setdiff(selected, numbers)
  if (length(unknown) > 0) stop("No such step(s): ", paste(unknown, collapse = ", "), "\n", usage, call. = FALSE)
}

# ---- Settings and modes --------------------------------------------------------------
# Sourcing run_config.R runs its checks, so a setting that cannot work stops the
# run here and not in the middle of a step.
source(file.path(root, "run_config.R"))
modes_text <- paste(config$modes, collapse = ", ")
cat("Analysis modes of this run (config$modes in run_config.R): ", modes_text, "\n", sep = "")

# Steps that repeat the analysis of the draft, or compare the draft with it, and
# need the legacy mode. Steps 09, 11 and 13 check this themselves when they are
# started alone and end without doing anything. Step 10 has no such check: run
# without the legacy mode it finishes and reports every term of the draft as
# carried by no eligible gene, which is wrong, so it is not started here either.
needs_legacy <- c(9L, 10L, 11L, 13L)
no_legacy <- !"legacy" %in% config$modes

# ---- Run ---------------------------------------------------------------------------------
run_steps <- function() {
  old_wd <- setwd(root); on.exit(setwd(old_wd), add = TRUE)
  rscript <- file.path(R.home("bin"), "Rscript")
  timing <- data.frame(step = character(0), minutes = numeric(0), result = character(0))
  for (i in which(numbers %in% selected)) {
    if (no_legacy && numbers[i] %in% needs_legacy) {
      cat("\n", basename(scripts[i]), " is not run: it needs the legacy mode, and \"legacy\" is not in config$modes (",
          modes_text, ").\n", sep = "")
      timing <- rbind(timing, data.frame(step = basename(scripts[i]), minutes = 0, result = "not run: needs the legacy mode"))
      next
    }
    t0 <- Sys.time()
    status <- system2(rscript, shQuote(scripts[i]))
    minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    timing <- rbind(timing, data.frame(step = basename(scripts[i]), minutes = round(minutes, 1),
                                       result = if (status == 0) "finished" else paste("failed, exit status", status)))
    if (status != 0) {
      print(timing)
      stop("Step failed: ", basename(scripts[i]), " (exit status ", status, ")", call. = FALSE)
    }
  }
  timing
}
timing <- run_steps()
not_run <- timing$step[timing$result != "finished"]
cat("\n", if (length(not_run) == 0) "All requested steps finished." else
      paste0("The requested steps finished, except those that need the legacy mode and were not run: ",
             paste(not_run, collapse = ", "), "."),
    " Analysis modes: ", modes_text, ".\n", sep = "")
print(timing)
