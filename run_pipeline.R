# run_pipeline.R
#
# Orchestrator. Runs the numbered scripts in 01_pipeline/ in order, each in a
# fresh R process, and stops at the first failure.
#
#   Rscript run_pipeline.R                        the paper: every step of the paper run, nothing else
#   Rscript run_pipeline.R --steps 5:7            steps 05, 06, 07 of the run
#   Rscript run_pipeline.R --steps 2:3,6:10       steps 02, 03 and 06 to 10 of the run
#   Rscript run_pipeline.R --list                 show the steps and exit
#   Rscript run_pipeline.R --plan [switches]      print what a run would do, and exit
#
# Optional analyses (NOT ADOPTED in the paper; each off by default). Switches
# (a) to (g) together define ONE analysis, written under
# 04_output/not_adopted/<id>/, where <id> joins the tokens of the switches that
# are on in the order corrected, cov5, mnchisq, ctrl or inter, ctexcl, d50, q05:
#   --coverage-filter        (a) token corrected: the coverage filter of the Methods
#                                (99.9th percentile of coverage per library; every line read)
#   --coverage 5             (b) token cov5: 5 reads per CpG in each library (--coverage 10 is the paper's value,
#                                and does not turn off coverage_standard = 5L set in run_config.R)
#   --overdispersion         (c) token mnchisq: methylKit overdispersion = "MN", test = "Chisq"
#   --exclude-ct-variants    (d) token ctexcl: CpGs at C-to-T variants removed (needs the list of
#                                01_pipeline/optional/ct_variant_loci.R; the run refuses without it)
#   --difference 50          (e) token d50: DML at a difference of more than 50 points instead of 20, as in
#                                Venkataraman et al. 2022 (--difference 20 is the paper's value, and does not turn
#                                off dml_difference_standard = 50L set in run_config.R); the fit of the analysis
#                                without d50 is reused
#   --qvalue 0.05            (f) token q05: DML at a q-value below 0.05 instead of 0.01 (--qvalue 0.01 is the
#                                paper's value, and does not turn off dml_qvalue_standard = 0.05 set in
#                                run_config.R); the fit of the analysis without q05 is reused
#   --controls-only-ploidy   (g) token ctrl: the ploidy contrast fitted on the 12 control libraries only (diploid
#                                against triploid controls, no covariate; issue I-07); the CpGs, the normalisation
#                                and the heat contrast are those of the analysis without ctrl (its fit is reused)
#   --interaction            (g) token inter: a tested ploidy by heat interaction (indicator triploid and
#                                heat-treated, covariates ploidy and heat; issue I-02) and its DML list, next to
#                                the ploidy and heat contrasts of the analysis without inter (its fit is reused).
#                                --controls-only-ploidy and --interaction cannot be given together (setting
#                                contrast_variant of run_config.R)
#   --gene-level-go          (h) step 10 and the gene-level part of step 11, for the analysis of the run
#   --diagnostics [all|name,...]  (i) scripts of 06_diagnostics/ (listed in 06_diagnostics/diagnostics.csv); they
#                                run after step 11 and before step 13 whatever --steps selects
# Examples: --coverage-filter -> not_adopted/corrected/; --coverage 5 --overdispersion ->
# not_adopted/cov5-mnchisq/; --overdispersion --difference 50 -> not_adopted/mnchisq-d50/;
# --qvalue 0.05 --difference 50 -> not_adopted/d50-q05/; --controls-only-ploidy -> not_adopted/ctrl/;
# --coverage-filter --coverage 5 --overdispersion --interaction -> not_adopted/corrected-cov5-mnchisq-inter/;
# --gene-level-go alone -> step 10 on the paper's DML, in not_adopted/legacy/.
#
# Run kinds. A run with no optional setting on (in run_config.R or by a switch)
# is the paper run: steps 01 to 09, 11 (paper part), 12 and 13 (paper part),
# writing the paper's outputs (04_output/<step>/, 05_manuscript/). Any switch makes
# an optional run, which never writes into a paper folder: steps 04 to 09 for an
# optional analysis, 10 and 11 (gene-level part) with (h), the diagnostics with
# (i), and step 13 (optional part, which reports values of analyses not run as
# "not run"). Steps 01 to 03 and 12 belong to the paper run only. The settings
# of run_config.R (list `optional`) and the switches are combined: a switch can
# only turn a setting on.
#
# The selection after --steps is one word without spaces. Anything else on the
# command line stops the run before a step starts. The run refuses to start when
# WGBS_ANALYSIS, WGBS_GENE_LEVEL_GO or WGBS_DIAGNOSTICS is set in its
# environment: they choose the analysis of a step started alone, and this script
# sets them for its steps (and unsets them when it ends, also after a failure).
#
# Settings are in run_config.R. They are read and checked here before the first
# step starts, and the plan of the run is printed at the start and the run kind
# and analysis at the end. The first run downloads 10.5 GB from gannet (step 01;
# the inputs of steps 04 and 05) and the 31 MB GO.db archive of step 09. Step 04
# needs about 9.5 GB of memory and step 05 about 10 GB with WGBS_CORES=4.
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

# Steps of each run kind. The paper run never starts step 10; an optional run
# never starts steps 01 to 03 and 12.
paper_steps    <- c(1:9, 11:13)
analysis_steps <- 4:9           # an optional analysis other than the paper's
gene_steps     <- 10:11         # with --gene-level-go
variant_script <- file.path(root, "01_pipeline", "optional", "ct_variant_loci.R")
diag_table     <- file.path(root, "06_diagnostics", "diagnostics.csv")

# ---- Command line ------------------------------------------------------------------
usage <- paste0(
  "Usage: Rscript run_pipeline.R [--steps <selection>] [--coverage-filter] [--coverage 5|10] [--overdispersion]\n",
  "                              [--controls-only-ploidy | --interaction]\n",
  "                              [--exclude-ct-variants] [--difference 20|50] [--qvalue 0.01|0.05]\n",
  "                              [--gene-level-go] [--diagnostics [all|name,...]] [--plan] [--list] [--help]\n",
  "  <selection> is one word without spaces: step numbers and ranges joined by commas, for example 5:7 or 3,12 or 2:3,6:10.\n",
  "  Without a switch the run is the paper's. Each switch turns on an optional analysis that the paper did not adopt;\n",
  "  its outputs go to 04_output/not_adopted/ and 05_manuscript/not_adopted/. --coverage 10, --difference 20 and\n",
  "  --qvalue 0.01 name the paper's values and do not turn off a setting made in run_config.R. See the header of run_pipeline.R.")
fail <- function(...) stop(..., "\n", usage, call. = FALSE)
opt <- list(steps = NULL, coverage_filter = FALSE, coverage = NULL, overdispersion = FALSE, controls_only = FALSE, interaction = FALSE,
            exclude_ct = FALSE, difference = NULL, qvalue = NULL, gene_level = FALSE, diagnostics = NULL, plan = FALSE, list = FALSE, help = FALSE)
seen <- character(0)
i <- 1
while (i <= length(args)) {
  a <- args[i]
  if (a %in% seen) fail("The switch ", a, " is given twice.")
  seen <- c(seen, a)
  takes_value <- function() { if (i == length(args) || startsWith(args[i + 1], "--")) fail("The switch ", a, " needs a value."); args[i + 1] }
  if (a %in% c("--help", "-h")) { opt$help <- TRUE; i <- i + 1 }
  else if (a == "--list") { opt$list <- TRUE; i <- i + 1 }
  else if (a == "--plan") { opt$plan <- TRUE; i <- i + 1 }
  else if (a == "--coverage-filter") { opt$coverage_filter <- TRUE; i <- i + 1 }
  else if (a == "--overdispersion") { opt$overdispersion <- TRUE; i <- i + 1 }
  else if (a == "--controls-only-ploidy") { opt$controls_only <- TRUE; i <- i + 1 }
  else if (a == "--interaction") { opt$interaction <- TRUE; i <- i + 1 }
  else if (a == "--exclude-ct-variants") { opt$exclude_ct <- TRUE; i <- i + 1 }
  else if (a == "--gene-level-go") { opt$gene_level <- TRUE; i <- i + 1 }
  else if (a == "--coverage") {
    v <- takes_value()
    if (!v %in% c("5", "10")) fail("--coverage must be followed by 5 or 10 (the paper); it is followed by '", v, "'.")
    opt$coverage <- as.integer(v); i <- i + 2
  } else if (a == "--difference") {
    v <- takes_value()
    if (!v %in% c("20", "50")) fail("--difference must be followed by 50 or 20 (the paper); it is followed by '", v, "'.")
    opt$difference <- as.integer(v); i <- i + 2
  } else if (a == "--qvalue") {
    v <- takes_value()
    if (!v %in% c("0.01", "0.05")) fail("--qvalue must be followed by 0.05 or 0.01 (the paper); it is followed by '", v, "'.")
    opt$qvalue <- v; i <- i + 2
  } else if (a == "--diagnostics") {
    if (i < length(args) && !startsWith(args[i + 1], "--")) { opt$diagnostics <- strsplit(args[i + 1], ",", fixed = TRUE)[[1]]; i <- i + 2 }
    else { opt$diagnostics <- "all"; i <- i + 1 }
  } else if (a == "--steps") {
    v <- takes_value()
    part <- "[0-9]{1,2}(:[0-9]{1,2})?"
    if (!grepl(sprintf("^%s(,%s)*$", part, part), v)) fail("Cannot read the selection after --steps: ", shQuote(v))
    opt$steps <- unlist(lapply(strsplit(v, ",")[[1]], function(p) {
      r <- as.integer(strsplit(p, ":")[[1]])
      if (length(r) == 2 && r[1] > r[2]) fail("A range must run upward: ", p)
      if (length(r) == 2) seq(r[1], r[2]) else r
    }))
    unknown <- setdiff(opt$steps, numbers)
    if (length(unknown) > 0) fail("No such step(s): ", paste(unknown, collapse = ", "))
    i <- i + 2
  } else fail("Cannot read the argument(s): ", paste(shQuote(args), collapse = " "))
}
if (opt$help) { cat(usage, "\n", sep = ""); quit(status = 0) }
if (opt$controls_only && opt$interaction) {
  fail("--controls-only-ploidy and --interaction cannot be given together: an analysis has the control-only ploidy contrast ",
       "(token ctrl) or the tested interaction (token inter), not both.")
}
if (opt$list) {
  kinds <- vapply(numbers, function(n) paste(c(if (n %in% paper_steps) "paper", if (n %in% analysis_steps) "optional analysis",
                                                  if (n %in% gene_steps) "--gene-level-go", if (n == 13) "optional (values not run)"),
                                                collapse = ", "), "")
  cat(sprintf("%02d  %-34s %s\n", numbers, basename(scripts), kinds), sep = "")
  cat("--  optional/ct_variant_loci.R       the C-to-T variant list for --exclude-ct-variants; run by hand (01_pipeline/optional/README.md)\n")
  if (file.exists(diag_table)) {
    d <- utils::read.csv(diag_table, stringsAsFactors = FALSE)
    cat(sprintf("--  06_diagnostics/%-28s --diagnostics%s\n", paste0(d$script, ".R"), ifelse(d$needs_network, " (needs the network)", "")), sep = "")
  }
  quit(status = 0)
}

# ---- Settings, switches and the plan -------------------------------------------------------
# Sourcing run_config.R runs its checks, so a setting that cannot work stops the
# run here and not in the middle of a step.
env_vars <- c("WGBS_ANALYSIS", "WGBS_GENE_LEVEL_GO", "WGBS_DIAGNOSTICS")
set_vars <- env_vars[!is.na(Sys.getenv(env_vars, unset = NA))]
if (length(set_vars) > 0) {
  stop("The environment variable(s) ", paste(set_vars, collapse = ", "), " are set in this shell. They choose the analysis of a ",
       "step started alone, and a variable left over would change the run without saying so. Unset them ",
       "(unset ", paste(set_vars, collapse = " "), ") and choose the analysis with run_config.R or a switch.", call. = FALSE)
}
old_wd <- setwd(root)
source(file.path(root, "run_config.R"))
source(file.path(root, "02_R_functions", "analyses.R"))   # analysis_id(), ct_variant_list_present()
o <- config$optional
tokens <- c(if (o$coverage_filter_of_methods || opt$coverage_filter) "corrected",
            if (o$coverage_standard == 5 || identical(opt$coverage, 5L)) "cov5",
            if (o$overdispersion_chisq || opt$overdispersion) "mnchisq",
            if (o$contrast_variant == "ctrl" || opt$controls_only) "ctrl",
            if (o$contrast_variant == "inter" || opt$interaction) "inter",
            if (o$exclude_ct_variants || opt$exclude_ct) "ctexcl",
            if (o$dml_difference_standard == 50 || identical(opt$difference, 50L)) "d50",
            if (o$dml_qvalue_standard == 0.05 || identical(opt$qvalue, "0.05")) "q05")
gene_level <- o$gene_level_go || opt$gene_level
diagnostics <- unique(c(o$diagnostics, opt$diagnostics))
if ("all" %in% diagnostics) diagnostics <- "all"
Sys.setenv(WGBS_ANALYSIS = analysis_id(tokens), WGBS_GENE_LEVEL_GO = as.character(gene_level),
           WGBS_DIAGNOSTICS = if (length(diagnostics)) paste(diagnostics, collapse = ",") else "none")
# The three variables are unset however the rest of this script ends (a failed
# step, a refused setting), so that a run sourced in an R session leaves none set.
# tryCatch(finally) and not on.exit(): at the top level of a sourced file
# on.exit() runs at once, and under Rscript it never runs. The block is
# evaluated in this environment, so its assignments stay visible.
tryCatch({
source(file.path(root, "run_config.R"))          # checks the combination; sets config$analysis, $addons, $run_kind
setwd(old_wd)
a <- config$analysis; kind <- config$run_kind

steps_of_kind <- if (kind == "paper") paper_steps else
  sort(c(if (!a$paper) analysis_steps, if (gene_level) gene_steps, 13L))
selected <- if (is.null(opt$steps)) steps_of_kind else opt$steps
diag_list <- character(0)
if (length(diagnostics) > 0) {
  d <- utils::read.csv(diag_table, stringsAsFactors = FALSE); d <- d[order(d$order), ]
  diag_list <- if (identical(diagnostics, "all")) d$script else d$script[d$script %in% diagnostics]
}
folder <- if (kind == "paper") "04_output/<step>/ and 05_manuscript/ (the paper)" else
  sprintf("04_output/not_adopted/%s/, 04_output/not_adopted/comparisons/ and 05_manuscript/not_adopted/", a$id)
variant_list <- ct_variant_list_present()
in_cache <- function(fit_id) file.exists(file.path(root, "04_output", "05_differential_methylation", "cache", sprintf("methylkit_%s.rds", fit_id)))
fit_cached <- in_cache(a$fit_id)
# Options ctrl and inter: step 05 takes the united table and the contrasts it
# does not change from the fit of the analysis without the token (base_fit_id),
# and fits only the new contrast, by sequence.
fit_line <- function() {
  base_name <- if (a$base_fit_id == "legacy") "the paper's fit" else sprintf("the fit of the analysis %s", a$base_fit_id)
  if (a$contrast_variant == "none") {
    return(sprintf("Fit of step 05: methylkit_%s.rds%s, %s", a$fit_id,
                   if (a$fit_id == "legacy") " (the paper's fit)" else if (a$fit_id != a$id) sprintf(" (the fit of the analysis %s)", a$fit_id) else "",
                   if (fit_cached) "in the cache (reused if its fingerprint is unchanged)" else "not in the cache (fitted by step 05)"))
  }
  new_contrast <- if (a$contrast_variant == "ctrl") "the control-only ploidy contrast" else "the interaction contrast"
  reused <- if (a$contrast_variant == "ctrl") "the united table and the heat contrast" else "the united table and the ploidy and heat contrasts"
  minutes <- if (a$min_coverage == 5) "about 45 minutes (an estimate)" else "about 12 minutes (an estimate)"
  whole <- if (a$min_coverage == 5) "about 90 minutes, by sequence (measured for the base fits)" else "about 25 minutes (measured for the base fits)"
  status <- if (fit_cached) {
    sprintf("in the cache (reused if its fingerprint, which holds that of methylkit_%s.rds, is unchanged)", a$base_fit_id)
  } else if (in_cache(a$base_fit_id)) {
    sprintf("not in the cache: fitted by step 05, %s by sequence, %s at WGBS_CORES=2 (methylkit_%s.rds is in the cache)",
            new_contrast, minutes, a$base_fit_id)
  } else {
    sprintf("not in the cache, nor is methylkit_%s.rds: step 05 fits %s first (%s at WGBS_CORES=2), then %s by sequence (%s)",
            a$base_fit_id, a$base_fit_id, whole, new_contrast, minutes)
  }
  sprintf("Fit of step 05: methylkit_%s.rds (%s; %s from methylkit_%s.rds, %s), %s", a$fit_id, new_contrast, reused,
          a$base_fit_id, base_name, status)
}
reason_not <- function(n) {
  if (kind == "paper" && n == 10) return("step 10 is not part of the paper; use --gene-level-go")
  if (kind == "optional" && n %in% c(1:3, 12)) return("not part of an optional run (it does not depend on the analysis)")
  if (kind == "optional" && a$paper && n %in% analysis_steps) return("the analysis is the paper's: its outputs are the paper's")
  if (kind == "optional" && n %in% gene_steps && !gene_level) return("needs --gene-level-go")
  "not selected"
}

plan_lines <- c(
  sprintf("Run kind: %s", if (kind == "paper") "paper (the results of the paper and nothing else)" else "OPTIONAL (not adopted in the paper)"),
  sprintf("Analysis: %s%s", a$id, if (a$paper) " (the paper's)" else
    sprintf(" (filter %s, %gx in every library, overdispersion %s%s, contrasts %s, C-to-T variants %s, DML at a difference of more than %g points and a q-value below %g)",
            a$filter, a$min_coverage, a$overdispersion,
            if (a$overdispersion == "MN") " with the chi-squared test" else "",
            switch(a$contrast_variant, none = "of the paper", ctrl = "with the ploidy contrast on the controls only (ctrl)",
                   inter = "of the paper and a tested interaction (inter)"),
            if (a$ct_excluded) "excluded" else "kept",
            a$dml_difference, a$dml_qvalue)),
  if (!a$paper) sprintf("Contrasts: %s", contrast_text(a)),
  sprintf("Outputs: %s", folder),
  sprintf("Steps to run: %s", paste(sprintf("%02d", intersect(selected, steps_of_kind)), collapse = " ")),
  if (length(diag_list)) sprintf("Diagnostics, after step 11 and before step 13, whatever --steps selects: %s", paste(diag_list, collapse = ", ")),
  if (length(setdiff(numbers, intersect(selected, steps_of_kind))))
    paste0("Not run: ", paste(sprintf("%02d (%s)", setdiff(numbers, intersect(selected, steps_of_kind)),
                                      vapply(setdiff(numbers, intersect(selected, steps_of_kind)), reason_not, "")), collapse = "; ")),
  sprintf("C-to-T variant list (04_output/not_adopted/ct_variants/tables/): %s", if (variant_list) "present" else "absent"),
  if (!a$paper) fit_line(),
  if (!a$paper) paste("Time and memory: at 10x (measured 2026-11-20) step 04 7 to 9 minutes, 9.2 GB; step 05 22 to 24 minutes per new fit at",
                       "WGBS_CORES=2 (about 10 GB with its workers; under a minute from the cache); steps 06 to 09 about 4 minutes. At 5x",
                       "(measured 2026-11-21) step 04 8 to 15 minutes; step 05 88 to 107 minutes per new fit at WGBS_CORES=2, fitted by sequence",
                       "(peak 9.3 GB; about a minute from the cache); steps 06 to 09 about 5 minutes"),
  if ("cov5" %in% a$tokens && !fit_cached && (a$contrast_variant == "none" || !in_cache(a$base_fit_id))) paste("NOTE: an analysis at 5x tests about 2.3 million CpGs instead of 0.48 million. Step 05 fits it by sequence",
                                  "and saves each part: about 3.2 hours at WGBS_CORES=1 and 1.6 hours at WGBS_CORES=2 (measured on one",
                                  "linkage group on 2026-11-20; peak 8.3 GB in the import, 6.6 GB in a test on 2 cores). A run that is",
                                  "stopped goes on from the saved parts when the same command is run again."),
  if ("cov5" %in% a$tokens && !fit_cached && a$contrast_variant != "none" && in_cache(a$base_fit_id))
    paste("NOTE: an analysis at 5x tests about 2.3 million CpGs instead of 0.48 million. Step 05 fits the new contrast by sequence",
          "on the united table of the cached fit and saves each part (about 45 minutes at WGBS_CORES=2, an estimate: one of the two",
          "contrasts of a fit). A run that is stopped goes on from the saved parts when the same command is run again."),
  if ("cov5" %in% a$tokens && fit_cached) paste("NOTE: an analysis at 5x tests about 2.3 million CpGs instead of 0.48 million. Its fit is in the cache",
                                 "and is reused (no refit), but steps 04 and 06 to 09 still handle about 2.3 million CpGs, and step 05",
                                 "reads the fit of all of them and writes their per-CpG tables (up to 98 MB each)."),
  if (gene_level) "Step 10 about 3 minutes; step 11 (gene-level part) about 2 minutes")
cat("Plan of this run\n", paste0("  ", plan_lines, collapse = "\n"), "\n", sep = "")
if (opt$plan) quit(status = 0)

if (a$ct_excluded && !variant_list) {
  stop("--exclude-ct-variants needs the C-to-T variant list, which is absent (04_output/not_adopted/ct_variants/tables/). Make it first with:\n  ",
       ct_variant_build_command, call. = FALSE)
}

# ---- Run ---------------------------------------------------------------------------------
diagnostic_ready <- function(script) {
  d <- utils::read.csv(diag_table, stringsAsFactors = FALSE)
  if (isTRUE(d$needs_network[d$script == script])) {
    ok <- tryCatch({ curlGetHeaders("https://gannet.fish.washington.edu/"); TRUE }, error = function(e) FALSE)
    if (!ok) return("needs the network (gannet.fish.washington.edu)")
  }
  # The outputs of the analysis 'corrected' that the diagnostic reads (column
  # reads_corrected_files), not only the stamp of its folder.
  if (grepl("corrected", d$reads_analyses[d$script == script])) {
    folder <- file.path(root, "04_output", "not_adopted", "corrected")
    needed <- c("analysis.csv", strsplit(d$reads_corrected_files[d$script == script], ";", fixed = TRUE)[[1]])
    absent <- needed[!file.exists(file.path(folder, needed))]
    if (length(absent) > 0) {
      return(sprintf("needs the outputs of the analysis 'corrected' (--coverage-filter); absent from 04_output/not_adopted/corrected/: %s",
                     paste(absent, collapse = ", ")))
    }
  }
  ""
}
run_steps <- function() {
  old_wd <- setwd(root); on.exit(setwd(old_wd), add = TRUE)
  rscript <- file.path(R.home("bin"), "Rscript")
  timing <- data.frame(step = character(0), minutes = numeric(0), result = character(0))
  run_one <- function(path, label) {
    t0 <- Sys.time()
    status <- system2(rscript, shQuote(path))
    minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    data.frame(step = label, minutes = round(minutes, 1), result = if (status == 0) "finished" else paste("failed, exit status", status))
  }
  for (n in sort(unique(c(numbers, if (length(diag_list)) 11.5)))) {
    if (n == 11.5) {                                  # the diagnostics, after step 11 and before step 13
      for (s in diag_list) {
        why <- diagnostic_ready(s)
        row <- if (nzchar(why)) data.frame(step = paste0("06_diagnostics/", s, ".R"), minutes = 0, result = paste("not run:", why))
               else run_one(file.path(root, "06_diagnostics", paste0(s, ".R")), paste0("06_diagnostics/", s, ".R"))
        timing <- rbind(timing, row)
        if (startsWith(row$result, "failed")) { print(timing); stop("Diagnostic failed: ", s, call. = FALSE) }
      }
      next
    }
    script <- scripts[numbers == n]
    if (!n %in% selected && !n %in% steps_of_kind) next
    if (!n %in% steps_of_kind || !n %in% selected) {
      timing <- rbind(timing, data.frame(step = basename(script), minutes = 0, result = paste("not run:", reason_not(n))))
      next
    }
    row <- run_one(script, basename(script))
    timing <- rbind(timing, row)
    if (row$result != "finished") {
      print(timing)
      stop("Step failed: ", basename(script), " (", row$result, ")", call. = FALSE)
    }
  }
  timing
}
timing <- run_steps()
not_run <- timing$step[timing$result != "finished"]
cat("\n", if (length(not_run) == 0) "All requested steps finished." else
      paste0("The requested steps finished; not run: ", paste(not_run, collapse = ", "), "."),
    sprintf(" Run kind: %s; analysis: %s%s.\n", kind, a$id,
            if (kind == "optional") sprintf(" (outputs under 04_output/not_adopted/%s/)", a$id) else " (the paper)"), sep = "")
print(timing)
}, finally = Sys.unsetenv(env_vars))
