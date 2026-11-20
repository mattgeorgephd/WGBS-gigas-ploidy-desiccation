# tests/check_fit_by_sequence.R
#
# The fit by sequence of step 05 (every optional analysis) against the fit in
# one piece (fit_mode(), the paper's fit). It writes into the cache and the
# output folders, so it runs only in a scratch copy of the repository without
# .git (tests/scratch_guard.R says how to make one), which must hold the cached
# fits in one piece of the analyses it compares
# (04_output/05_differential_methylation/cache/methylkit_<id>.rds; by default
# legacy and mnchisq). From the root of the copy:
#   CT_TEST_SCRATCH=1 WGBS_CACHE_DIR=<download cache> WGBS_CORES=4 Rscript tests/check_fit_by_sequence.R [id ...]
# For each analysis (about 15 minutes each at 10x on 4 cores):
#  1. the cached fit in one piece is current: its fingerprint is fit_fingerprint(a) now;
#  2. fit_by_sequence() (import and every part fitted and saved) returns an object
#     identical() to it: filtering record, united table, and both methylDiff
#     objects (p-values, q-values, differences, row names, slots);
#  3. write_mode() of that fit writes every table and BED file of the analysis
#     with the content it had before (the paper's tracked files for legacy; for
#     an optional analysis the files of its folder, the per-CpG tables included
#     when they are on disk); .gz files are compared decompressed;
#  4. resume: a second call fits nothing (every saved part is reused) and gives
#     the identical object; after one part is deleted, one is cut short and a
#     file of another fingerprint is added, a third call fits those two parts
#     again and nothing else, deletes the third file, and gives the identical
#     object;
#  5. an edit of a function of the fit by sequence changes fit_fingerprint(a, TRUE)
#     and leaves fit_fingerprint(a) (the fit in one piece) as it was.
# Exit status 1 when a check fails. The parts are left in the copy's cache.

suppressPackageStartupMessages(library(here))
source(here("tests", "scratch_guard.R"))
scratch_guard("tests/check_fit_by_sequence.R")
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "methylKit"))

# The definitions of step 05 (its functions, `step` and `lib`), without running the step.
for (e in parse(here("01_pipeline", "05_differential_methylation.R"), keep.source = FALSE)) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      ((is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) || as.character(e[[2]]) %in% c("step", "lib"))) {
    eval(e, globalenv())
  }
}

ids <- commandArgs(trailingOnly = TRUE)
if (length(ids) == 0) ids <- c("legacy", "mnchisq")
fails <- character(0)
check <- function(ok, what) {
  cat(if (isTRUE(ok)) "ok  " else "FAIL", what, "\n")
  if (!isTRUE(ok)) fails <<- c(fails, what)
  invisible(isTRUE(ok))
}
md5_content <- function(f) {
  if (!endsWith(f, ".gz")) return(unname(tools::md5sum(f)))
  tmp <- tempfile(); on.exit(unlink(tmp))
  con <- gzfile(f, "rb"); out <- file(tmp, "wb")
  repeat { b <- readBin(con, "raw", 1e7); if (length(b) == 0) break; writeBin(b, out) }
  close(con); close(out)
  unname(tools::md5sum(tmp))
}
# The filtering record is a data.table, whose attribute .internal.selfref is a
# pointer that differs between a table in memory and one read from a file (where
# it is nil); it is compared after a round trip through serialize().
round_trip <- function(o) unserialize(serialize(o, NULL))
same_fit <- function(x, y) {
  c(mode = identical(x$mode, y$mode), filtering = identical(round_trip(x$filtering), round_trip(y$filtering)),
    united = identical(x$united, y$united),
    ploidy = identical(x$diffmeth$ploidy, y$diffmeth$ploidy), heat = identical(x$diffmeth$heat, y$diffmeth$heat))
}
describe_difference <- function(x, y) {
  for (contrast in c("ploidy", "heat")) {
    a <- getData(x$diffmeth[[contrast]]); b <- getData(y$diffmeth[[contrast]])
    if (nrow(a) != nrow(b)) { cat("    ", contrast, ": rows", nrow(a), "and", nrow(b), "\n"); next }
    for (col in c("pvalue", "qvalue", "meth.diff")) {
      cat(sprintf("     %s %s: %d of %d values differ, largest absolute difference %g\n", contrast, col,
                  sum(a[[col]] != b[[col]]), nrow(a), max(abs(a[[col]] - b[[col]]))))
    }
  }
}

for (id in ids) {
  a <- analysis_spec(id)
  cat("\n== analysis", id, "(fit id", a$fit_id, ")\n")
  config$analysis <- a
  config$run_kind <- if (a$paper) "paper" else "optional"
  whole_file <- file.path(step_cache_dir(step), sprintf("methylkit_%s.rds", a$fit_id))
  if (!check(file.exists(whole_file), paste("1. cached fit in one piece present:", basename(whole_file)))) next
  whole <- readRDS(whole_file)
  check(identical(whole$fingerprint, fit_fingerprint(a)), "1. it is current (fingerprint of the fit in one piece unchanged)")

  fp <- fit_fingerprint(a, TRUE)
  unlink(parts_dir(a), recursive = TRUE)
  t1 <- Sys.time()
  fit <- fit_by_sequence(a, fp)
  cat(sprintf("     fit by sequence: %.1f min on %d core(s)\n", as.numeric(difftime(Sys.time(), t1, units = "mins")), config$cores))
  same <- same_fit(fit, whole)
  check(all(same), paste("2. identical to the fit in one piece:", paste(names(same), same, sep = " ", collapse = ", ")))
  if (!all(same)) describe_difference(fit, whole)
  for (contrast in c("ploidy", "heat")) {
    d <- getData(fit$diffmeth[[contrast]])
    cat(sprintf("     %s: %d CpGs, %d p-values of exactly 1, none NA: %s; DML %d\n", contrast, nrow(d), sum(d$pvalue == 1),
                !anyNA(d$pvalue), sum(d$qvalue < config$dml_qvalue & abs(d$meth.diff) > config$dml_difference)))
  }

  # 3. The written files.
  folders <- c(analysis_dir(a, step, "tables"), analysis_dir(a, step, "bed"))
  files <- list.files(folders, full.names = TRUE)
  before <- vapply(files, md5_content, character(1))
  invisible(write_mode(c(fit, list(fingerprint = fp)), a))
  after_files <- list.files(folders, full.names = TRUE)
  after <- vapply(files, md5_content, character(1))
  changed <- basename(files)[before != after]
  check(length(changed) == 0 && setequal(after_files, files),
        sprintf("3. write_mode() rewrote %d files with the same content%s%s", length(files),
                if (length(changed)) paste0("; changed: ", paste(changed, collapse = ", ")) else "",
                if (!setequal(after_files, files)) paste0("; new: ", paste(basename(setdiff(after_files, files)), collapse = ", ")) else ""))

  # 4. Resume.
  part_files <- function() list.files(parts_dir(a), full.names = TRUE)
  mtimes <- function() setNames(file.mtime(part_files()), basename(part_files()))
  m0 <- mtimes()
  Sys.sleep(1.1)
  fit2 <- fit_by_sequence(a, fp)
  m1 <- mtimes()
  check(identical(m0, m1) && all(same_fit(fit2, fit)), sprintf("4. second call: all %d saved files reused, nothing fitted, identical result", length(m0)))
  parts <- grep("_(ploidy|heat)_", names(m0), value = TRUE)
  deleted <- parts[2]; cut <- parts[length(parts)]
  unlink(file.path(parts_dir(a), deleted))
  writeBin(readBin(file.path(parts_dir(a), cut), "raw", 100), file.path(parts_dir(a), cut))
  stranger <- file.path(parts_dir(a), paste0("000000000000_ploidy_", sub("^[0-9a-f]{12}_ploidy_", "", parts[1])))
  file.copy(file.path(parts_dir(a), parts[1]), stranger)
  Sys.sleep(1.1)
  m1 <- mtimes()
  fit3 <- fit_by_sequence(a, fp)
  m2 <- mtimes()
  refitted <- names(m2)[is.na(m1[names(m2)]) | m2 != m1[names(m2)]]
  check(setequal(refitted, c(deleted, cut)) && !file.exists(stranger) && setequal(names(m2), names(m0)) && all(same_fit(fit3, fit)),
        sprintf("4. third call: refitted %s (expected %s and %s), file of another fingerprint deleted: %s, identical result",
                paste(refitted, collapse = " and "), deleted, cut, !file.exists(stranger)))
  whole_fp <- whole$fingerprint
  rm(whole, fit, fit2, fit3); invisible(gc(verbose = FALSE))

  # 5. Fingerprints.
  keep <- test_part
  test_part <- function(obj, contrast, a) { x <- keep(obj, contrast, a); x }
  check(!identical(fit_fingerprint(a, TRUE), fp) && identical(fit_fingerprint(a), whole_fp),
        "5. an edit of test_part() changes the fingerprint of the fit by sequence, not that of the fit in one piece")
  test_part <- keep
  check(identical(fit_fingerprint(a, TRUE), fp), "5. with test_part() restored, the fingerprint is the one of the fit")
}

cat("\n", if (length(fails)) paste(length(fails), "check(s) failed") else "all checks pass", "\n")
quit(status = if (length(fails)) 1 else 0)
