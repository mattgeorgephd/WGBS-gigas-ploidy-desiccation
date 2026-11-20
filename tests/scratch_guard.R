# tests/scratch_guard.R
#
# Sourced first by the two tests that delete or overwrite files of the
# repository (tests/check_step13_not_run.R, tests/check_ct_exclusion.R). They
# must run in a scratch copy WITHOUT .git, never in the original checkout:
#   mkdir <scratch> && tar --exclude=./.git --exclude=./04_output/not_adopted/ct_variants -cf - . | tar -xf - -C <scratch>
#     (from the repository root; or cp -a, then delete <scratch>/.git and <scratch>/04_output/not_adopted/ct_variants/)
#   cd <scratch> && CT_TEST_SCRATCH=1 Rscript tests/<test>.R
# The guard stops the test, before it touches a file, when any of these holds:
#   - the working tree has a .git (folder or file), or git finds a remote for it
#     (it lies inside a clone);
#   - 04_output/not_adopted/ct_variants/tables/ct_variant_run.csv exists and is
#     not the synthetic file that check_ct_exclusion.R writes: the real C-to-T
#     variant list is there, and check_ct_exclusion.R would overwrite it;
#   - CT_TEST_SCRATCH is not set to 1.

scratch_guard <- function(test) {
  root <- here::here()
  why <- character(0)
  if (file.exists(file.path(root, ".git"))) why <- c(why, "the working tree has a .git")
  remotes <- suppressWarnings(tryCatch(system2("git", c("-C", shQuote(root), "remote"), stdout = TRUE, stderr = FALSE),
                                       error = function(e) character(0)))
  if (length(remotes) > 0) why <- c(why, paste0("git finds a remote for it (", paste(remotes, collapse = ", "), ")"))
  run_file <- file.path(root, "04_output", "not_adopted", "ct_variants", "tables", "ct_variant_run.csv")
  if (file.exists(run_file) && !synthetic_run_file(run_file)) {
    why <- c(why, "04_output/not_adopted/ct_variants/tables/ct_variant_run.csv is the record of the real variant screen")
  }
  if (!identical(Sys.getenv("CT_TEST_SCRATCH"), "1")) why <- c(why, "CT_TEST_SCRATCH=1 is not set")
  if (length(why) > 0) {
    stop(test, " deletes and overwrites files and runs only in a scratch copy of the repository without .git and without ",
         "04_output/not_adopted/ct_variants/ (see tests/scratch_guard.R). Refused: ", paste(why, collapse = "; "), ".", call. = FALSE)
  }
  invisible(TRUE)
}
# The run record that check_ct_exclusion.R writes next to its synthetic list.
synthetic_run_marker <- data.frame(section = "synthetic", library_id = "", item = "written_by", value = "tests/check_ct_exclusion.R")
synthetic_run_file <- function(f) {
  r <- tryCatch(utils::read.csv(f, colClasses = "character", na.strings = NULL), error = function(e) NULL)
  is.data.frame(r) && identical(dim(r), c(1L, 4L)) && identical(unname(unlist(r[1, ])), unname(unlist(synthetic_run_marker[1, ])))
}
