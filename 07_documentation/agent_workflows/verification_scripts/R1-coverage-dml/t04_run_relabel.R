# Run the real step 04 script, expression by expression, on regional extracts, with all paths redirected to scratch.
suppressPackageStartupMessages(library(here))
.h <- new.env()
.h$repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
.h$scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/R1-coverage-dml"
stopifnot(normalizePath(here()) == .h$repo)
.h$heads_dir <- Sys.getenv("T04_HEADS", file.path(.h$scratch, "t05/heads")); .h$out_root <- Sys.getenv("T04_OUT", file.path(.h$scratch, "t04/out_relabel"))
unlink(.h$out_root, recursive = TRUE)
.h$run <- function() {
  ex <- parse(file.path(.h$repo, "01_pipeline/04_cpg_coverage.R"))
  for (k in seq_along(ex)) {
    eval(ex[[k]], globalenv())
    if (any(grepl("02_R_functions", deparse(ex[[k]])))) {
      assign("fetch_gannet_file", function(file, config, ...) if (grepl("\\.fa$", file)) file.path(.h$scratch, "t04/test_genome.fa") else file.path(.h$heads_dir, file), globalenv())
      assign("output_dir", function(step, sub = NULL) { d <- if (is.null(sub)) file.path(.h$out_root, step) else file.path(.h$out_root, step, sub); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }, globalenv())
      assign("output_path", function(step, sub, file) file.path(output_dir(step, sub), file), globalenv())
      assign("step_cache_dir", function(step) output_dir(step, "cache"), globalenv())
      .h$orig_read <- read_libraries; assign("read_libraries", function() { lib <- .h$orig_read(); lib$treatment[lib$zymo_project == "zr3534"] <- NA; lib }, globalenv()); cat("paths redirected; sample sheet as read_libraries() returns it when heat rows are labelled \"heat treated\"\n")
    }
  }
}
.h$run()
save(cov_mat, nC_mat, universe, chr_levels, first_row, cutoffs, file = file.path(.h$out_root, "t04_objects.RData"))
