# run_pipeline.R
#
# Orchestrator. Runs the numbered scripts in 01_pipeline/ in order, each in a
# fresh R process, and stops at the first failure.
#
#   Rscript run_pipeline.R                 all steps
#   Rscript run_pipeline.R --steps 5:7     steps 05, 06, 07
#   Rscript run_pipeline.R --steps 3,12    steps 03 and 12
#   Rscript run_pipeline.R --list          show the steps and exit
#
# Settings are in run_config.R. Steps 04 and 05 download 9.8 GB from gannet the
# first time (via step 01) and need about 12 GB of memory.

suppressPackageStartupMessages(library(here))
args <- commandArgs(trailingOnly = TRUE)

scripts <- sort(list.files(here("01_pipeline"), pattern = "^[0-9]{2}_.*\\.R$", full.names = TRUE))
numbers <- as.integer(substr(basename(scripts), 1, 2))

if ("--list" %in% args) {
  cat(sprintf("%02d  %s\n", numbers, basename(scripts)), sep = "")
  quit(status = 0)
}
selected <- numbers
if ("--steps" %in% args) {
  spec <- args[which(args == "--steps") + 1]
  selected <- unlist(lapply(strsplit(spec, ",")[[1]], function(p) {
    r <- as.integer(strsplit(p, ":")[[1]]); if (length(r) == 2) seq(r[1], r[2]) else r
  }))
  unknown <- setdiff(selected, numbers)
  if (length(unknown) > 0) stop("No such step(s): ", paste(unknown, collapse = ", "))
} else if (length(args) > 0) {
  stop("Unknown argument(s): ", paste(args, collapse = " "), ". Use --steps or --list.")
}

rscript <- file.path(R.home("bin"), "Rscript")
timing <- data.frame(step = character(0), minutes = numeric(0))
for (i in which(numbers %in% selected)) {
  t0 <- Sys.time()
  status <- system2(rscript, shQuote(scripts[i]))
  minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  timing <- rbind(timing, data.frame(step = basename(scripts[i]), minutes = round(minutes, 1)))
  if (status != 0) {
    print(timing)
    stop("Step failed: ", basename(scripts[i]), " (exit status ", status, ")")
  }
}
cat("\nAll requested steps finished.\n"); print(timing)
