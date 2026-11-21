# tests/check_cached_fit.R
#
# The cache logic of step 05 (cached_fit(), the lock of a fit by sequence, the
# removal of leftovers, and the key of the saved parts), with stand-ins for the
# fitting functions, so it fits nothing and takes about a minute. It moves and
# writes files of the cache folder (and restores the fits it moves), so it runs
# only in a scratch copy of the repository without .git (tests/scratch_guard.R)
# that holds the cached fits methylkit_legacy.rds and methylkit_mnchisq.rds of
# 04_output/05_differential_methylation/cache/. From the root of the copy:
#   CT_TEST_SCRATCH=1 Rscript tests/check_cached_fit.R
# Checks:
#  1. legacy: the cached fit is used (its fingerprint is unchanged), no fitting
#     function is called, and leftovers are removed: parts_legacy/ and a
#     temporary file of a process that has ended; a temporary file of a live
#     process is kept;
#  2. legacy: a cached fit whose fingerprint has the two extra fields of a fit by
#     sequence is not used; the fit in one piece is made again;
#  3. mnchisq: a cached fit made in one piece is used by the fit by sequence
#     (fallback), and its leftover parts_mnchisq/ is removed;
#  4. ctexcl (fit id legacy) uses methylkit_legacy.rds; without it, ctexcl is
#     fitted in one piece, never by sequence;
#  5. lock: a lock of a process that has ended, and one whose process id now
#     belongs to another process (other start time), are taken over; while the
#     fit runs the lock names this process; a lock of a live process stops the
#     step with a message that names it; the lock is gone afterwards;
#  6. the key of the saved parts is the same under LC_ALL=C and LC_ALL=C.UTF-8
#     (the md5 of digest_md5(), used for the stamp only, is not).
# Exit status 1 when a check fails.

suppressPackageStartupMessages(library(here))
args <- commandArgs(trailingOnly = TRUE)
print_key <- length(args) == 2 && args[1] == "--print-key"
if (!print_key) {
  source(here("tests", "scratch_guard.R"))
  scratch_guard("tests/check_cached_fit.R")
}
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "methylKit"))
for (e in parse(here("01_pipeline", "05_differential_methylation.R"), keep.source = FALSE)) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      ((is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) || as.character(e[[2]]) %in% c("step", "lib"))) {
    eval(e, globalenv())
  }
}
if (print_key) {       # used by check 6, in a child process with another locale
  fp <- fit_fingerprint(analysis_spec(args[2]), TRUE)
  cat(Sys.getlocale("LC_CTYPE"), fingerprint_key(fp), digest_md5(fp), "\n")
  quit(status = 0)
}

fails <- character(0)
check <- function(ok, what) {
  cat(if (isTRUE(ok)) "ok  " else "FAIL", what, "\n")
  if (!isTRUE(ok)) fails <<- c(fails, what)
}
cache <- step_cache_dir(step)
cache_file <- function(id) file.path(cache, sprintf("methylkit_%s.rds", id))
calls <- character(0)
stub <- function(kind) function(a, ...) {
  calls <<- c(calls, kind)
  if (kind == "parts") {
    lock <- file.path(parts_dir(a), "lock")
    calls <<- c(calls, if (file.exists(lock) && identical(read_lock(lock)$pid, Sys.getpid())) "locked by me" else "not locked")
  }
  list(mode = a$fit_id, filtering = data.frame(), united = NULL, diffmeth = list(), stand_in = kind)
}
run <- function(id) {
  calls <<- character(0)
  tryCatch(cached_fit(analysis_spec(id), fit_whole = stub("whole"), fit_parts = stub("parts")),
           error = function(e) structure(conditionMessage(e), class = "failed"))
}
dead_pid <- function() {                          # the id of a process that has ended
  pid <- as.integer(system("sh -c 'echo $$'", intern = TRUE))
  stopifnot(!isTRUE(tools::pskill(pid, 0L)))
  pid
}
for (id in c("legacy", "mnchisq")) stopifnot(file.exists(cache_file(id)))

# 1. legacy, valid, with leftovers.
dir.create(parts_dir(analysis_spec("legacy")), showWarnings = FALSE)
writeLines("x", file.path(parts_dir(analysis_spec("legacy")), "0123456789ab_ploidy_X.rds"))
gone <- sprintf("%s.%d.part", cache_file("legacy"), dead_pid()); writeLines("x", gone)
kept <- sprintf("%s.%d.part", cache_file("legacy"), Sys.getpid()); writeLines("x", kept)
fit <- run("legacy")
check(!inherits(fit, "failed") && is.null(fit$stand_in) && length(calls) == 0,
      "1. legacy: the cached fit is used, nothing is fitted (the paper's fingerprint is unchanged)")
check(!dir.exists(parts_dir(analysis_spec("legacy"))) && !file.exists(gone) && file.exists(kept),
      "1. leftovers removed: parts_legacy/ and the temporary file of an ended process; that of a live process kept")
unlink(kept)
paper_fp <- fit$fingerprint
rm(fit); invisible(gc(verbose = FALSE))

# 2. legacy with the extra fields of a fit by sequence.
backup <- paste0(cache_file("legacy"), ".backup")
stopifnot(file.rename(cache_file("legacy"), backup))
tryCatch({
  x <- readRDS(backup)
  x$fingerprint <- fit_fingerprint(analysis_spec("legacy"), TRUE)
  saveRDS(x, cache_file("legacy")); rm(x)
  fit <- run("legacy")
  check(identical(fit$stand_in, "whole") && identical(calls, "whole"),
        "2. legacy: a cached fit with the fingerprint of a fit by sequence is not used; the fit in one piece is made again")
}, finally = { unlink(cache_file("legacy")); file.rename(backup, cache_file("legacy")) })

# 3. mnchisq: fallback to its fit in one piece; leftovers removed.
dir.create(parts_dir(analysis_spec("mnchisq")), showWarnings = FALSE)
writeLines("x", file.path(parts_dir(analysis_spec("mnchisq")), "0123456789ab_united.rds"))
fit <- run("mnchisq")
check(!inherits(fit, "failed") && is.null(fit$stand_in) && length(calls) == 0 &&
        identical(fit$fingerprint, fit_fingerprint(analysis_spec("mnchisq"), FALSE)),
      "3. mnchisq: the cached fit made in one piece is used by the fit by sequence; nothing is fitted")
check(!dir.exists(parts_dir(analysis_spec("mnchisq"))), "3. leftover parts_mnchisq/ removed")
rm(fit); invisible(gc(verbose = FALSE))

# 4. ctexcl.
fit <- run("ctexcl")
check(!inherits(fit, "failed") && is.null(fit$stand_in) && length(calls) == 0 && identical(fit$fingerprint, paper_fp),
      "4. ctexcl uses methylkit_legacy.rds (the paper's fit)")
rm(fit); invisible(gc(verbose = FALSE))
stopifnot(file.rename(cache_file("legacy"), backup))
tryCatch({
  fit <- run("ctexcl")
  check(identical(calls, "whole") && !dir.exists(parts_dir(analysis_spec("ctexcl"))),
        "4. without methylkit_legacy.rds, ctexcl is fitted in one piece, not by sequence")
}, finally = { unlink(cache_file("legacy")); file.rename(backup, cache_file("legacy")) })

# 5. Lock, on a fit id that has no cached fit in the copy.
id <- "cov5-mnchisq"
a <- analysis_spec(id)
stopifnot(!file.exists(cache_file(id)))
lock <- file.path(parts_dir(a), "lock")
dir.create(parts_dir(a), recursive = TRUE, showWarnings = FALSE)
writeLines(paste(dead_pid(), "12345", "2026-01-01 00:00:00 UTC"), lock)
fit <- run(id)
check(identical(calls, c("parts", "locked by me")) && !file.exists(lock) && !dir.exists(parts_dir(a)),
      "5. a lock of an ended process is taken over; the lock names this process during the fit and is gone afterwards")
unlink(cache_file(id))
sleeper <- as.integer(system("sh -c 'sleep 120 > /dev/null 2>&1 & echo $!'", intern = TRUE))
dir.create(parts_dir(a), recursive = TRUE, showWarnings = FALSE)
writeLines(paste(sleeper, "1", "2026-01-01 00:00:00 UTC"), lock)
fit <- run(id)
check(identical(calls, c("parts", "locked by me")), "5. a lock whose process id now belongs to another process (other start time) is taken over")
unlink(cache_file(id))
dir.create(parts_dir(a), recursive = TRUE, showWarnings = FALSE)
writeLines(paste(sleeper, process_start(sleeper), "2026-01-01 00:00:00 UTC"), lock)
fit <- run(id)
check(inherits(fit, "failed") && grepl("Another run of step 05", fit) && grepl(as.character(sleeper), fit) && length(calls) == 0 &&
        file.exists(lock) && !file.exists(cache_file(id)),
      paste0("5. a lock of a live process stops the step, nothing fitted, lock left: \"", substr(fit, 1, 120), "...\""))
tools::pskill(sleeper)
unlink(parts_dir(a), recursive = TRUE)

# 6. The key under two locales.
out <- vapply(c("C", "C.UTF-8"), function(loc) {
  system2("env", c(paste0("LC_ALL=", loc), "Rscript", "tests/check_cached_fit.R", "--print-key", "cov5"), stdout = TRUE, stderr = FALSE) |>
    tail(1)
}, character(1))
parts <- strsplit(trimws(out), " ")
cat("     ", out, sep = "\n      ")
check(length(parts) == 2 && all(lengths(parts) == 3) && parts[[1]][2] == parts[[2]][2] &&
        parts[[1]][2] == fingerprint_key(fit_fingerprint(analysis_spec("cov5"), TRUE)),
      "6. the key of the saved parts is the same under LC_ALL=C and LC_ALL=C.UTF-8, and in this session")
cat(sprintf("      (digest_md5, the stamp's record only: %s under C and %s under C.UTF-8)\n",
            if (length(parts) == 2) parts[[1]][3] else "?", if (length(parts) == 2) parts[[2]][3] else "?"))

cat("\n", if (length(fails)) paste(length(fails), "check(s) failed") else "all checks pass", "\n")
quit(status = if (length(fails)) 1 else 0)
