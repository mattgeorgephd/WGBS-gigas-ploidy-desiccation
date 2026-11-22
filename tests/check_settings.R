# tests/check_settings.R
#
# T0 of the settings layout (07_documentation/agent_workflows/scope_inputs/settings_build_spec.md,
# section 8): static and unit tests of the settings, the switches and the write
# guards. Seconds; reads no data. Run from the repository root:
#   Rscript tests/check_settings.R [BASE commit, default 66c7121]
# Prints one line per check and ends with exit status 1 if any check fails.

suppressPackageStartupMessages({ library(here); library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
base <- if (length(args)) args[1] else "66c7121"
results <- list()
ok <- function(name, pass, detail = "") {
  results[[length(results) + 1]] <<- data.table(check = name, pass = isTRUE(pass), detail = detail)
  cat(sprintf("%-4s %s%s\n", if (isTRUE(pass)) "ok" else "FAIL", name, if (nzchar(detail)) paste0("  (", detail, ")") else ""))
}
env_vars <- c("WGBS_ANALYSIS", "WGBS_GENE_LEVEL_GO", "WGBS_DIAGNOSTICS")
Sys.unsetenv(env_vars)

# Settings as run_config.R gives them under the environment `env` (named values).
settings <- function(env = character(0), file = here("run_config.R")) {
  Sys.unsetenv(env_vars)
  if (length(env)) do.call(Sys.setenv, as.list(env))
  on.exit(Sys.unsetenv(env_vars))
  e <- new.env()
  source(file, local = e)
  e$config
}
stops <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")

# ---- Static ---------------------------------------------------------------------------
code <- c(list.files(c("01_pipeline", "02_R_functions", "06_diagnostics"), pattern = "\\.R$", full.names = TRUE),
          "run_config.R", "run_pipeline.R")
code <- code[!grepl("01_pipeline/optional/", code)]    # the variant step: not edited while its screen runs
hits <- unlist(lapply(code, function(f) { l <- readLines(f, warn = FALSE); f[grepl("config\\$modes|require_mode", l)] }))
ok("config$modes and require_mode() appear in no code file", length(hits) == 0, paste(unique(hits), collapse = ", "))
# Steps 05 to 11 apply the q cutoff of the analysis (a$dml_qvalue), never the frozen paper value.
steps_q <- list.files("01_pipeline", pattern = "^(0[5-9]|1[01])_.*\\.R$", full.names = TRUE)
hits_q <- unlist(lapply(steps_q, function(f) { l <- sub("#.*$", "", readLines(f, warn = FALSE)); f[grepl("config\\$dml_qvalue", l)] }))
ok("steps 05 to 11 use no config$dml_qvalue outside comments", length(steps_q) == 7 && length(hits_q) == 0, paste(unique(hits_q), collapse = ", "))
# Options ctrl and inter: the code that every existing fit depends on (the functions in
# the fingerprint of fit_fingerprint(), and fit_fingerprint() itself) is that of commit
# 76a900b (before the options), as R parses it, so the paper's fit and every cached fit
# keep their fingerprints.
fit_functions <- function(rev) {
  e <- new.env()
  files <- c("01_pipeline/05_differential_methylation.R", list.files("02_R_functions", pattern = "\\.R$", full.names = TRUE))
  for (f in files) {
    text <- if (is.null(rev)) readLines(f, warn = FALSE) else
      suppressWarnings(system2("git", c("show", paste0(rev, ":", f)), stdout = TRUE, stderr = FALSE))
    if (length(text) == 0) next
    for (x in parse(text = text, keep.source = FALSE)) {
      if (is.call(x) && identical(x[[1]], as.name("<-")) && is.name(x[[2]]) && is.call(x[[3]]) && identical(x[[3]][[1]], as.name("function"))) {
        assign(as.character(x[[2]]), paste(deparse(x[[3]]), collapse = "\n"), envir = e)
      }
    }
  }
  e
}
fit_code <- c("fit_mode", "read_filtered_counts", "read_cov", "first_line_lost", "upper_coverage_cutoff", "cpg_key", "key_to_position",
              "scale_counts_like_methylkit", "make_methylraw", "fit_by_sequence", "import_united", "sequence_parts", "test_part",
              "diffmeth_from_parts", "fit_fingerprint", "fingerprint_key")
now_f <- fit_functions(NULL); old_f <- fit_functions("76a900b")
changed_f <- fit_code[!vapply(fit_code, function(n) exists(n, now_f) && exists(n, old_f) && identical(get(n, now_f), get(n, old_f)), logical(1))]
ok("the fitting code of every existing fit (fit_mode, fit_by_sequence, their helpers, fit_fingerprint) is that of 76a900b",
   length(changed_f) == 0, paste(changed_f, collapse = ", "))

# ---- 192 ids x gene-level GO x diagnostics ----------------------------------------------------
source(here("02_R_functions", "analyses.R"))
source(here("02_R_functions", "paths.R"))
# The write paths are resolved with dir.create() masked: the functions above are
# defined in this (global) environment and find this version first, so the test
# creates no folder. Removed after the loop.
dir.create <- function(...) invisible(TRUE)
paper_dirs <- here(c("04_output", "05_manuscript"))
steps_written <- c("04_cpg_coverage", "05_differential_methylation", "06_dml_figures", "07_genome_features", "08_gene_annotation",
                   "09_go_enrichment_legacy", "10_go_enrichment_gene_level", "11_semantic_similarity")
under <- function(p, d) startsWith(p, paste0(d, "/"))
tokens6 <- c("corrected", "cov5", "mnchisq", "ctexcl", "d50", "q05")
ids6 <- unique(vapply(0:63, function(k) analysis_id(tokens6[bitwAnd(k, 2^(0:5)) > 0]), ""))
ids <- unique(c(ids6, unlist(lapply(c("ctrl", "inter"), function(v) vapply(0:63, function(k) analysis_id(c(tokens6[bitwAnd(k, 2^(0:5)) > 0], v)), "")))))
ok("192 analysis ids (64 without ctrl or inter, 64 with ctrl, 64 with inter), all valid",
   length(ids6) == 64 && length(ids) == 192 && identical(analysis_tokens, c("corrected", "cov5", "mnchisq", "ctrl", "inter", "ctexcl", "d50", "q05")) &&
     sum(grepl("(^|-)ctrl(-|$)", ids)) == 64 && sum(grepl("(^|-)inter(-|$)", ids)) == 64 &&
     !any(grepl("ctrl", ids) & grepl("inter", ids)) && all(vapply(ids, function(id) !stops(analysis_spec(id)), logical(1))),
   paste(length(ids), "ids"))
bad <- character(0)
for (id in ids) for (g in c("TRUE", "FALSE")) for (d in c("none", "all")) {
  cfg <- settings(c(WGBS_ANALYSIS = id, WGBS_GENE_LEVEL_GO = g, WGBS_DIAGNOSTICS = d))
  a <- cfg$analysis
  want_kind <- if (id == "legacy" && g == "FALSE" && d == "none") "paper" else "optional"
  fit_want <- analysis_id(setdiff(a$tokens, c("ctexcl", "d50", "q05")))
  variant_want <- if ("ctrl" %in% a$tokens) "ctrl" else if ("inter" %in% a$tokens) "inter" else "none"
  if (!(cfg$min_coverage == 10 && a$min_coverage == (if ("cov5" %in% a$tokens) 5 else 10) && a$fit_id == fit_want &&
        a$contrast_variant == variant_want && a$base_fit_id == analysis_id(setdiff(a$tokens, c("ctrl", "inter", "ctexcl", "d50", "q05"))) &&
        cfg$dml_difference == 20 && a$dml_difference == (if ("d50" %in% a$tokens) 50 else 20) &&
        cfg$dml_qvalue == 0.01 && a$dml_qvalue == (if ("q05" %in% a$tokens) 0.05 else 0.01) &&
        cfg$run_kind == want_kind && identical(cfg$addons$gene_level_go, g == "TRUE") &&
        identical(length(cfg$addons$diagnostics) > 0, d == "all"))) bad <- c(bad, paste(id, g, d))
  # paths of the analysis: never in a paper folder for an optional analysis, never under not_adopted/ for the paper in a paper run
  config <- cfg
  # (read path of a table the paper has: for the paper the paper folder, for another analysis its own folder)
  p <- tryCatch(analysis_input(a, "05_differential_methylation", "tables", "dml_counts_by_mode.csv"), error = function(e) "error")
  if (!a$paper && !under(p, not_adopted_path(a$id))) bad <- c(bad, paste("path", id))
  if (a$paper && p != here("04_output", "05_differential_methylation", "tables", "dml_counts_by_mode.csv")) bad <- c(bad, paste("paper path", id))
  # write paths: in an optional run every analysis (the paper's included) writes under 04_output/not_adopted/<id>/
  # and 05_manuscript/not_adopted/ only; in the paper run the paper writes its paper folders
  own <- not_adopted_path(a$id)
  for (st in steps_written) for (sub in c("tables", "figures")) {
    w <- c(tryCatch(analysis_dir(a, st, sub), error = function(e) "error"), tryCatch(analysis_path(a, st, sub, "x.csv"), error = function(e) "error"))
    if (cfg$run_kind == "optional" && !all(under(w, own))) bad <- c(bad, paste("write path", id, g, d, st))
    if (cfg$run_kind == "paper" && !all(w == c(here("04_output", st, sub), here("04_output", st, sub, "x.csv")))) bad <- c(bad, paste("paper write path", id, st))
  }
  m <- tryCatch(manuscript_path("not_adopted", a$id, "figures", "x.png"), error = function(e) "error")
  if (cfg$run_kind == "optional" && (!under(m, here("05_manuscript", "not_adopted")) || !stops(manuscript_path("figures", "x.png")))) {
    bad <- c(bad, paste("manuscript path", id, g, d))
  }
  if (cfg$run_kind == "optional" && !a$paper && stops(output_dir("05_differential_methylation", "tables")) == FALSE) bad <- c(bad, paste("output_dir", id, g, d))
}
rm(dir.create)
ok("768 combinations: min_coverage 10, dml_difference 20, dml_qvalue 0.01, analysis min_coverage, dml_difference, dml_qvalue, contrast_variant, fit_id, base_fit_id, run kind, add-ons, read and write paths",
   length(bad) == 0, paste(bad, collapse = "; "))
config <- settings()                                    # the paper run
ok("paper run: manuscript_path() without arguments is 05_manuscript/", identical(manuscript_path(), here("05_manuscript")))
ok("manuscript_path() rejects a '..' component", stops(manuscript_path("not_adopted", "..", "x.csv")) && stops(manuscript_path("figures/../../x.png")))
ok("fit_id of ctexcl is legacy; of cov5-mnchisq-ctexcl is cov5-mnchisq",
   analysis_spec("ctexcl")$fit_id == "legacy" && analysis_spec("cov5-mnchisq-ctexcl")$fit_id == "cov5-mnchisq")
ok("fit_id of d50 is legacy; of mnchisq-d50 mnchisq; of corrected-cov5-ctexcl-d50 corrected-cov5; of corrected-cov5-mnchisq-ctexcl-d50 corrected-cov5-mnchisq",
   analysis_spec("d50")$fit_id == "legacy" && analysis_spec("mnchisq-d50")$fit_id == "mnchisq" &&
     analysis_spec("corrected-cov5-ctexcl-d50")$fit_id == "corrected-cov5" &&
     analysis_spec("corrected-cov5-mnchisq-ctexcl-d50")$fit_id == "corrected-cov5-mnchisq")
ok("d50: dml_difference 50, every other setting that of the analysis without d50",
   all(vapply(setdiff(ids, ids[grepl("d50", ids)]), function(id) {
     x <- analysis_spec(id); y <- analysis_spec(analysis_id(c(x$tokens, "d50")))
     x$dml_difference == 20 && y$dml_difference == 50 &&
       identical(x[c("filter", "min_coverage", "overdispersion", "test", "ct_excluded", "dml_qvalue", "fit_id")],
                 y[c("filter", "min_coverage", "overdispersion", "test", "ct_excluded", "dml_qvalue", "fit_id")])
   }, logical(1))))
ok("fit_id of q05 is legacy; of mnchisq-q05 mnchisq; of d50-q05 legacy; of corrected-cov5-mnchisq-ctexcl-d50-q05 corrected-cov5-mnchisq",
   analysis_spec("q05")$fit_id == "legacy" && analysis_spec("mnchisq-q05")$fit_id == "mnchisq" &&
     analysis_spec("d50-q05")$fit_id == "legacy" &&
     analysis_spec("corrected-cov5-mnchisq-ctexcl-d50-q05")$fit_id == "corrected-cov5-mnchisq")
ok("q05: dml_qvalue 0.05, every other setting that of the analysis without q05",
   all(vapply(setdiff(ids, ids[grepl("q05", ids)]), function(id) {
     x <- analysis_spec(id); y <- analysis_spec(analysis_id(c(x$tokens, "q05")))
     x$dml_qvalue == 0.01 && y$dml_qvalue == 0.05 &&
       identical(x[c("filter", "min_coverage", "overdispersion", "test", "ct_excluded", "dml_difference", "fit_id")],
                 y[c("filter", "min_coverage", "overdispersion", "test", "ct_excluded", "dml_difference", "fit_id")])
   }, logical(1))))
ok("q05 is the last token of an id: corrected-cov5-mnchisq-ctexcl-d50-q05, mnchisq-q05, d50-q05",
   analysis_id(c("q05", "d50", "ctexcl", "mnchisq", "cov5", "corrected")) == "corrected-cov5-mnchisq-ctexcl-d50-q05" &&
     analysis_id(c("q05", "mnchisq")) == "mnchisq-q05" && analysis_id(c("q05", "d50")) == "d50-q05")
ok("ctrl and inter come after mnchisq and before ctexcl: corrected-cov5-mnchisq-ctrl-ctexcl-d50-q05, mnchisq-inter-q05, inter-ctexcl",
   analysis_id(c("q05", "d50", "ctexcl", "ctrl", "mnchisq", "cov5", "corrected")) == "corrected-cov5-mnchisq-ctrl-ctexcl-d50-q05" &&
     analysis_id(c("q05", "inter", "mnchisq")) == "mnchisq-inter-q05" && analysis_id(c("ctexcl", "inter")) == "inter-ctexcl")
ok("fit ids of ctrl and inter keep the token: ctrl -> ctrl, corrected-cov5-mnchisq-inter-ctexcl-d50-q05 -> corrected-cov5-mnchisq-inter, mnchisq-ctrl-q05 -> mnchisq-ctrl",
   analysis_spec("ctrl")$fit_id == "ctrl" && analysis_spec("corrected-cov5-mnchisq-inter-ctexcl-d50-q05")$fit_id == "corrected-cov5-mnchisq-inter" &&
     analysis_spec("mnchisq-ctrl-q05")$fit_id == "mnchisq-ctrl" && analysis_spec("inter-d50")$fit_id == "inter")
ok("base_fit_id: ctrl -> legacy, inter-ctexcl -> legacy, cov5-mnchisq-inter -> cov5-mnchisq, corrected-ctrl-q05 -> corrected; without ctrl or inter base_fit_id = fit_id",
   analysis_spec("ctrl")$base_fit_id == "legacy" && analysis_spec("inter-ctexcl")$base_fit_id == "legacy" &&
     analysis_spec("cov5-mnchisq-inter")$base_fit_id == "cov5-mnchisq" && analysis_spec("corrected-ctrl-q05")$base_fit_id == "corrected" &&
     all(vapply(ids6, function(id) { x <- analysis_spec(id); x$base_fit_id == x$fit_id && x$contrast_variant == "none" }, logical(1))))
ok("ctrl and inter: every other setting that of the analysis without the token (filter, coverage, model, ctexcl, thresholds)",
   all(vapply(ids6, function(id) {
     x <- analysis_spec(id)
     all(vapply(c("ctrl", "inter"), function(v) {
       y <- analysis_spec(analysis_id(c(x$tokens, v)))
       y$contrast_variant == v && y$base_fit_id == x$fit_id && !y$paper &&
         identical(x[c("filter", "min_coverage", "overdispersion", "test", "ct_excluded", "dml_difference", "dml_qvalue")],
                   y[c("filter", "min_coverage", "overdispersion", "test", "ct_excluded", "dml_difference", "dml_qvalue")])
     }, logical(1)))
   }, logical(1))))
ok("an id with both ctrl and inter stops, saying they are mutually exclusive: ctrl-inter, corrected-cov5-mnchisq-ctrl-inter-q05",
   stops(analysis_spec("ctrl-inter")) && stops(analysis_spec("corrected-cov5-mnchisq-ctrl-inter-q05")) &&
     analysis_id(c("inter", "ctrl")) == "ctrl-inter" &&
     grepl("mutually exclusive", tryCatch(analysis_spec("ctrl-inter"), error = conditionMessage)))
# Stamps: fixtures in a temporary folder (not_adopted_path() is pointed there for these checks).
local({
  real_path <- not_adopted_path
  tmp <- tempfile("stamps"); dir.create(tmp)
  on.exit({ not_adopted_path <<- real_path; unlink(tmp, recursive = TRUE) })
  not_adopted_path <<- function(...) file.path(tmp, ...)
  stamp <- function(id, dml_difference = NULL, dml_qvalue = NULL, contrast_variant = NULL) {
    a <- analysis_spec(id)
    row <- data.table(id = a$id, filter = a$filter, min_coverage = format(a$min_coverage), overdispersion = a$overdispersion,
                      test = a$test)
    if (!is.null(contrast_variant)) row[, contrast_variant := contrast_variant]
    row[, ct_excluded := as.character(a$ct_excluded)]
    if (!is.null(dml_difference)) row[, dml_difference := dml_difference]
    if (!is.null(dml_qvalue)) row[, dml_qvalue := dml_qvalue]
    row[, `:=`(fit_id = a$fit_id, fit_fingerprint_md5 = "", ct_variant_list_md5 = "", git_commit = "", written_utc = "", steps_done = "04")]
    dir.create(file.path(tmp, id), showWarnings = FALSE)
    fwrite(row, file.path(tmp, id, "analysis.csv"))
  }
  stamp("mnchisq")                                       # as written before the option d50 existed
  s_old <- read_stamp("mnchisq")
  ok("a stamp written before d50 (no column dml_difference) reads with dml_difference 20, after ct_excluded, and passes check_stamp()",
     identical(s_old$dml_difference, "20") && identical(match("dml_difference", names(s_old)), match("ct_excluded", names(s_old)) + 1L) &&
       !stops(check_stamp(analysis_spec("mnchisq"), needs = "04")))
  stamp("d50", "20")
  ok("check_stamp() stops for the analysis d50 when its stamp says dml_difference 20", stops(check_stamp(analysis_spec("d50"), needs = "04")))
  stamp("mnchisq", "50")
  ok("check_stamp() stops for the analysis mnchisq when its stamp says dml_difference 50", stops(check_stamp(analysis_spec("mnchisq"), needs = "04")))
  stamp("mnchisq-d50", "50")
  ok("check_stamp() passes for mnchisq-d50 with a stamp that says dml_difference 50", !stops(check_stamp(analysis_spec("mnchisq-d50"), needs = "04")))
  ok("a stamp written before q05 (no column dml_qvalue) reads with dml_qvalue 0.01, after dml_difference, and passes check_stamp()",
     with(list(s = read_stamp("mnchisq-d50")), identical(s$dml_qvalue, "0.01") &&
            identical(match("dml_qvalue", names(s)), match("dml_difference", names(s)) + 1L) &&
            identical(names(s)[seq_along(stamp_columns)], stamp_columns)) &&
       !stops(check_stamp(analysis_spec("mnchisq-d50"), needs = "04")))
  stamp("ctexcl")                                        # written before d50 and q05: neither column
  ok("a stamp without dml_difference and dml_qvalue reads with 20 and 0.01 in the order of stamp_columns",
     with(list(s = read_stamp("ctexcl")), identical(s$dml_difference, "20") && identical(s$dml_qvalue, "0.01") &&
            identical(names(s)[seq_along(stamp_columns)], stamp_columns)))
  stamp("q05", "20", "0.01")
  ok("check_stamp() stops for the analysis q05 when its stamp says dml_qvalue 0.01", stops(check_stamp(analysis_spec("q05"), needs = "04")))
  stamp("mnchisq", "20", "0.05")
  ok("check_stamp() stops for the analysis mnchisq when its stamp says dml_qvalue 0.05", stops(check_stamp(analysis_spec("mnchisq"), needs = "04")))
  stamp("mnchisq-d50-q05", "50", "0.05")
  ok("check_stamp() passes for mnchisq-d50-q05 with a stamp that says dml_difference 50 and dml_qvalue 0.05",
     !stops(check_stamp(analysis_spec("mnchisq-d50-q05"), needs = "04")))
  # write_stamp() and write_analyses_included() write the column dml_qvalue after dml_difference
  paper_config <- config
  config <<- list(run_kind = "optional")
  on.exit(config <<- paper_config, add = TRUE)
  write_table <<- function(x, path) data.table::fwrite(x, path)
  on.exit(rm(write_table, envir = globalenv()), add = TRUE)
  ok("a stamp written before ctrl and inter (no column contrast_variant) reads with contrast_variant none, after test, and passes check_stamp()",
     with(list(s = read_stamp("mnchisq-d50-q05")), identical(s$contrast_variant, "none") &&
            identical(match("contrast_variant", names(s)), match("test", names(s)) + 1L) &&
            identical(names(s)[seq_along(stamp_columns)], stamp_columns)) &&
       !stops(check_stamp(analysis_spec("mnchisq-d50-q05"), needs = "04")))
  stamp("mnchisq-ctrl", "20", "0.01")                    # no column contrast_variant: reads as none
  ok("check_stamp() stops for mnchisq-ctrl when its stamp has no contrast_variant (read as none)",
     stops(check_stamp(analysis_spec("mnchisq-ctrl"), needs = "04")))
  stamp("mnchisq-ctrl", "20", "0.01", "inter")
  ok("check_stamp() stops for mnchisq-ctrl when its stamp says contrast_variant inter", stops(check_stamp(analysis_spec("mnchisq-ctrl"), needs = "04")))
  stamp("mnchisq-ctrl", "20", "0.01", "ctrl")
  ok("check_stamp() passes for mnchisq-ctrl with a stamp that says contrast_variant ctrl", !stops(check_stamp(analysis_spec("mnchisq-ctrl"), needs = "04")))
  stamp("inter-q05", "20", "0.05", "none")
  ok("check_stamp() stops for inter-q05 when its stamp says contrast_variant none", stops(check_stamp(analysis_spec("inter-q05"), needs = "04")))
  write_stamp(analysis_spec("cov5-q05"), "04_cpg_coverage")
  w <- fread(file.path(tmp, "cov5-q05", "analysis.csv"), colClasses = "character")
  write_analyses_included()
  inc <- fread(file.path(tmp, "comparisons", "analyses_included.csv"), colClasses = "character")
  ok("write_stamp() writes dml_qvalue 0.05 for cov5-q05 after dml_difference; analyses_included.csv carries dml_qvalue (0.01 for the paper and old stamps)",
     identical(w$dml_qvalue, "0.05") && identical(names(w)[seq_along(stamp_columns)], stamp_columns) &&
       identical(names(inc)[seq_along(stamp_columns)], stamp_columns) &&
       identical(inc[id == "legacy", dml_qvalue], "0.01") && identical(inc[id == "cov5-q05", dml_qvalue], "0.05") &&
       identical(inc[id == "ctexcl", dml_qvalue], "0.01") && identical(inc[id == "mnchisq-d50-q05", dml_qvalue], "0.05"))
  write_stamp(analysis_spec("cov5-inter-d50"), "04_cpg_coverage")
  w <- fread(file.path(tmp, "cov5-inter-d50", "analysis.csv"), colClasses = "character")
  write_analyses_included()
  inc <- fread(file.path(tmp, "comparisons", "analyses_included.csv"), colClasses = "character")
  ok("write_stamp() writes contrast_variant inter and fit_id cov5-inter for cov5-inter-d50; analyses_included.csv carries contrast_variant (none for the paper and old stamps)",
     identical(w$contrast_variant, "inter") && identical(w$fit_id, "cov5-inter") && identical(names(w)[seq_along(stamp_columns)], stamp_columns) &&
       identical(inc[id == "legacy", contrast_variant], "none") && identical(inc[id == "cov5-inter-d50", contrast_variant], "inter") &&
       identical(inc[id == "ctexcl", contrast_variant], "none") && identical(inc[id == "mnchisq-ctrl", contrast_variant], "ctrl") &&
       identical(names(inc)[seq_along(stamp_columns)], stamp_columns))
})

# ---- Invalid values stop -------------------------------------------------------------------
for (v in c("cov5-corrected", "legacy-cov5", "Cov5", "5x", " ", "d50-ctexcl", "d50-mnchisq", "d20", "D50", "legacy-d50",
            "q05-d50", "q05-mnchisq", "q01", "Q05", "q5", "legacy-q05", "q05-q05",
            "ctrl-inter", "inter-ctrl", "mnchisq-ctrl-inter", "ctrl-mnchisq", "ctexcl-inter", "inter-cov5", "CTRL", "Inter", "ctrl-ctrl",
            "legacy-ctrl", "interaction", "ctrl0")) {
  ok(sprintf("WGBS_ANALYSIS='%s' stops", v), stops(settings(c(WGBS_ANALYSIS = v))))
}
edited <- function(from, to) {
  f <- tempfile(fileext = ".R"); l <- readLines(here("run_config.R"))
  i <- grep(from, l, fixed = TRUE); stopifnot(length(i) == 1); l[i] <- sub(from, to, l[i], fixed = TRUE); writeLines(l, f); f
}
ok("coverage_standard = 7 stops", stops(settings(file = edited("coverage_standard          = 10L", "coverage_standard          = 7L"))))
ok("dml_difference_standard = 30 stops", stops(settings(file = edited("dml_difference_standard    = 20L", "dml_difference_standard    = 30L"))))
ok("dml_difference_standard = 50 gives the analysis d50 at 50 points, config$dml_difference still 20",
   with(settings(file = edited("dml_difference_standard    = 20L", "dml_difference_standard    = 50L")),
        analysis$id == "d50" && analysis$dml_difference == 50 && dml_difference == 20 && run_kind == "optional"))
f_q1 <- edited("dml_qvalue_standard        = 0.01,", "dml_qvalue_standard        = 0.1,")
ok("dml_qvalue_standard = 0.1 stops", stops(settings(file = f_q1)))
ok("dml_qvalue_standard = 0.05 gives the analysis q05 at q 0.05, config$dml_qvalue still 0.01, fit legacy",
   with(settings(file = edited("dml_qvalue_standard        = 0.01,", "dml_qvalue_standard        = 0.05,")),
        analysis$id == "q05" && analysis$dml_qvalue == 0.05 && dml_qvalue == 0.01 && analysis$fit_id == "legacy" && run_kind == "optional"))
f_q2 <- edited("dml_qvalue                = 0.01,", "dml_qvalue                = 0.05,")
ok("dml_qvalue = 0.05 stops (frozen)", stops(settings(file = f_q2)))
cv <- 'contrast_variant           = "none",'
ok("contrast_variant = \"both\" stops", stops(settings(file = edited(cv, 'contrast_variant           = "both",'))))
ok("contrast_variant = NA stops", stops(settings(file = edited(cv, 'contrast_variant           = NA,'))))
ok("contrast_variant = c(\"ctrl\", \"inter\") stops", stops(settings(file = edited(cv, 'contrast_variant           = c("ctrl", "inter"),'))))
ok("contrast_variant = \"ctrl\" gives the analysis ctrl, fit ctrl on the base fit legacy",
   with(settings(file = edited(cv, 'contrast_variant           = "ctrl",')),
        analysis$id == "ctrl" && analysis$fit_id == "ctrl" && analysis$base_fit_id == "legacy" && run_kind == "optional"))
f_inter <- edited(cv, 'contrast_variant           = "inter",')
ok("contrast_variant = \"inter\" gives the analysis inter; with WGBS_ANALYSIS=ctrl the variable wins (ctrl)",
   with(settings(file = f_inter), analysis$id == "inter" && analysis$contrast_variant == "inter") &&
     settings(c(WGBS_ANALYSIS = "ctrl"), file = f_inter)$analysis$id == "ctrl")
ok("min_coverage = 5 stops (frozen)", stops(settings(file = edited("min_coverage              = 10,", "min_coverage              = 5,"))))
ok("dml_difference = 15 stops (frozen)", stops(settings(file = edited("dml_difference            = 20,", "dml_difference            = 15,"))))
ok("WGBS_GENE_LEVEL_GO=yes stops", stops(settings(c(WGBS_GENE_LEVEL_GO = "yes"))))
ok("WGBS_DIAGNOSTICS=nosuchscript stops", stops(settings(c(WGBS_DIAGNOSTICS = "nosuchscript"))))
rscript <- file.path(R.home("bin"), "Rscript")
out <- suppressWarnings(system2(rscript, c("run_pipeline.R", "--plan"), stdout = TRUE, stderr = TRUE, env = "WGBS_ANALYSIS=cov5"))
ok("run_pipeline.R refuses to start with WGBS_ANALYSIS set", !is.null(attr(out, "status")) && any(grepl("are set in this shell", out)))

# ---- Write guards ---------------------------------------------------------------------------
config <- settings(c(WGBS_ANALYSIS = "cov5"))
ok("optional run: output_dir(05, tables) stops", stops(output_dir("05_differential_methylation", "tables")))
ok("optional run: manuscript_path(figures) stops", stops(manuscript_path("figures")))
ok("optional run: output_dir(05, cache) allowed", !stops(output_dir("05_differential_methylation", "cache")))
ok("optional run: manuscript_path(not_adopted, ...) allowed", !stops(manuscript_path("not_adopted", "x.csv")))
config <- settings()
ok("paper run: analysis_dir of another analysis stops", stops(analysis_dir(analysis_spec("cov5"), "04_cpg_coverage", "tables")))
rm(config)

# ---- --plan of the examples ------------------------------------------------------------------
plans <- list(
  list(args = character(0), id = "legacy", steps = "01 02 03 04 05 06 07 08 09 11 12 13"),
  list(args = "--coverage-filter", id = "corrected", steps = "04 05 06 07 08 09 13"),
  list(args = c("--coverage", "5"), id = "cov5", steps = "04 05 06 07 08 09 13"),
  list(args = "--overdispersion", id = "mnchisq", steps = "04 05 06 07 08 09 13"),
  list(args = "--exclude-ct-variants", id = "ctexcl", steps = "04 05 06 07 08 09 13"),
  list(args = "--gene-level-go", id = "legacy", steps = "10 11 13"),
  list(args = "--diagnostics", id = "legacy", steps = "13"),
  list(args = c("--coverage-filter", "--coverage", "5", "--exclude-ct-variants"), id = "corrected-cov5-ctexcl", steps = "04 05 06 07 08 09 13"),
  list(args = c("--coverage", "5", "--overdispersion", "--exclude-ct-variants"), id = "cov5-mnchisq-ctexcl", steps = "04 05 06 07 08 09 13"),
  list(args = c("--coverage", "5", "--overdispersion"), id = "cov5-mnchisq", steps = "04 05 06 07 08 09 13"),
  list(args = c("--overdispersion", "--gene-level-go"), id = "mnchisq", steps = "04 05 06 07 08 09 10 11 13"),
  list(args = c("--difference", "20"), id = "legacy", steps = "01 02 03 04 05 06 07 08 09 11 12 13"),
  list(args = c("--difference", "50"), id = "d50", steps = "04 05 06 07 08 09 13"),
  list(args = c("--overdispersion", "--difference", "50"), id = "mnchisq-d50", steps = "04 05 06 07 08 09 13"),
  list(args = c("--difference", "50", "--coverage", "5", "--exclude-ct-variants"), id = "cov5-ctexcl-d50", steps = "04 05 06 07 08 09 13"),
  list(args = c("--coverage-filter", "--coverage", "5", "--overdispersion", "--exclude-ct-variants"), id = "corrected-cov5-mnchisq-ctexcl",
       steps = "04 05 06 07 08 09 13"),
  list(args = c("--coverage-filter", "--coverage", "5", "--overdispersion", "--exclude-ct-variants", "--difference", "50"),
       id = "corrected-cov5-mnchisq-ctexcl-d50", steps = "04 05 06 07 08 09 13"),
  list(args = c("--difference", "50", "--gene-level-go"), id = "d50", steps = "04 05 06 07 08 09 10 11 13"),
  list(args = c("--qvalue", "0.01"), id = "legacy", steps = "01 02 03 04 05 06 07 08 09 11 12 13"),
  list(args = c("--qvalue", "0.05"), id = "q05", steps = "04 05 06 07 08 09 13"),
  list(args = c("--qvalue", "0.05", "--difference", "50"), id = "d50-q05", steps = "04 05 06 07 08 09 13"),
  list(args = c("--overdispersion", "--qvalue", "0.05"), id = "mnchisq-q05", steps = "04 05 06 07 08 09 13"),
  list(args = c("--qvalue", "0.05", "--coverage", "5", "--exclude-ct-variants"), id = "cov5-ctexcl-q05", steps = "04 05 06 07 08 09 13"),
  list(args = c("--coverage-filter", "--coverage", "5", "--overdispersion", "--exclude-ct-variants", "--difference", "50", "--qvalue", "0.05"),
       id = "corrected-cov5-mnchisq-ctexcl-d50-q05", steps = "04 05 06 07 08 09 13"),
  list(args = c("--qvalue", "0.05", "--gene-level-go"), id = "q05", steps = "04 05 06 07 08 09 10 11 13"),
  list(args = "--controls-only-ploidy", id = "ctrl", steps = "04 05 06 07 08 09 13"),
  list(args = "--interaction", id = "inter", steps = "04 05 06 07 08 09 13"),
  list(args = c("--interaction", "--overdispersion"), id = "mnchisq-inter", steps = "04 05 06 07 08 09 13"),
  list(args = c("--qvalue", "0.05", "--controls-only-ploidy", "--coverage", "5"), id = "cov5-ctrl-q05", steps = "04 05 06 07 08 09 13"),
  list(args = c("--coverage-filter", "--coverage", "5", "--overdispersion", "--interaction", "--exclude-ct-variants", "--difference", "50", "--qvalue", "0.05"),
       id = "corrected-cov5-mnchisq-inter-ctexcl-d50-q05", steps = "04 05 06 07 08 09 13"),
  list(args = c("--controls-only-ploidy", "--gene-level-go"), id = "ctrl", steps = "04 05 06 07 08 09 10 11 13"))
for (p in plans) {
  out <- suppressWarnings(system2(rscript, c("run_pipeline.R", "--plan", p$args), stdout = TRUE, stderr = TRUE))
  got_id <- sub("^  Analysis: ([^ ]+).*", "\\1", grep("^  Analysis:", out, value = TRUE))
  got_steps <- sub("^  Steps to run: ", "", grep("^  Steps to run:", out, value = TRUE))
  ok(sprintf("--plan %s", paste(p$args, collapse = " ")), identical(got_id, p$id) && identical(got_steps, p$steps),
     paste(got_id, "|", got_steps))
}
out <- suppressWarnings(system2(rscript, c("run_pipeline.R", "--coverage", "7"), stdout = TRUE, stderr = TRUE))
ok("--coverage 7 stops before a step starts", !is.null(attr(out, "status")))
out <- suppressWarnings(system2(rscript, c("run_pipeline.R", "--difference", "30"), stdout = TRUE, stderr = TRUE))
ok("--difference 30 stops before a step starts", !is.null(attr(out, "status")))
out <- suppressWarnings(system2(rscript, c("run_pipeline.R", "--plan", "--overdispersion", "--difference", "50"), stdout = TRUE, stderr = TRUE))
ok("--plan of mnchisq-d50 names 50 points and the fit of mnchisq",
   any(grepl("DML at a difference of more than 50 points", out, fixed = TRUE)) && any(grepl("Fit of step 05: methylkit_mnchisq.rds", out, fixed = TRUE)))
out <- suppressWarnings(system2(rscript, c("run_pipeline.R", "--qvalue", "0.1", "--plan"), stdout = TRUE, stderr = TRUE))
ok("--qvalue 0.1 stops before a step starts", !is.null(attr(out, "status")))
out <- suppressWarnings(system2(rscript, c("run_pipeline.R", "--qvalue", ".05", "--plan"), stdout = TRUE, stderr = TRUE))
ok("--qvalue .05 stops before a step starts", !is.null(attr(out, "status")))
out <- suppressWarnings(system2(rscript, c("run_pipeline.R", "--plan", "--overdispersion", "--difference", "50", "--qvalue", "0.05"), stdout = TRUE, stderr = TRUE))
ok("--plan of mnchisq-d50-q05 names 50 points, q 0.05 and the fit of mnchisq",
   any(grepl("Analysis: mnchisq-d50-q05 (", out, fixed = TRUE)) &&
     any(grepl("DML at a difference of more than 50 points and a q-value below 0.05", out, fixed = TRUE)) &&
     any(grepl("Fit of step 05: methylkit_mnchisq.rds (the fit of the analysis mnchisq)", out, fixed = TRUE)))
out <- suppressWarnings(system2(rscript, c("run_pipeline.R", "--plan", "--qvalue", "0.05"), stdout = TRUE, stderr = TRUE))
ok("--plan of q05 names q 0.05 and the paper's fit",
   any(grepl("DML at a difference of more than 20 points and a q-value below 0.05", out, fixed = TRUE)) &&
     any(grepl("Fit of step 05: methylkit_legacy.rds (the paper's fit)", out, fixed = TRUE)))
out <- suppressWarnings(system2(rscript, c("run_pipeline.R", "--plan", "--qvalue", "0.01", "--overdispersion"), stdout = TRUE, stderr = TRUE))
ok("--plan --qvalue 0.01 --overdispersion is mnchisq at q 0.01",
   any(grepl("Analysis: mnchisq (", out, fixed = TRUE)) && any(grepl("a q-value below 0.01", out, fixed = TRUE)))
out <- suppressWarnings(system2(rscript, c("run_pipeline.R", "--overdispersion", "--overdispersion", "--plan"), stdout = TRUE, stderr = TRUE))
ok("a switch given twice stops", !is.null(attr(out, "status")))
out <- suppressWarnings(system2(rscript, c("run_pipeline.R", "--plan", "--controls-only-ploidy", "--interaction"), stdout = TRUE, stderr = TRUE))
ok("--controls-only-ploidy with --interaction stops before a step starts, saying they cannot be given together",
   !is.null(attr(out, "status")) && any(grepl("cannot be given together", out, fixed = TRUE)) && !any(grepl("Plan of this run", out, fixed = TRUE)))
out <- suppressWarnings(system2(rscript, c("run_pipeline.R", "--plan", "--controls-only-ploidy"), stdout = TRUE, stderr = TRUE))
ok("--plan of ctrl names the control-only contrast, the fit ctrl, the heat contrast of the paper's fit, and the state of the cache",
   any(grepl("Analysis: ctrl (", out, fixed = TRUE)) && any(grepl("contrasts with the ploidy contrast on the controls only (ctrl)", out, fixed = TRUE)) &&
     any(grepl("Contrasts: ploidy on the 12 control libraries only", out, fixed = TRUE)) &&
     any(grepl("Fit of step 05: methylkit_ctrl.rds (the control-only ploidy contrast; the united table and the heat contrast from methylkit_legacy.rds, the paper's fit), ", out, fixed = TRUE)) &&
     any(grepl("in the cache \\(reused|not in the cache: fitted by step 05, the control-only ploidy contrast by sequence, about 12 minutes|not in the cache, nor is methylkit_legacy.rds", out)))
out <- suppressWarnings(system2(rscript, c("run_pipeline.R", "--plan", "--interaction", "--coverage", "5", "--overdispersion", "--qvalue", "0.05"), stdout = TRUE, stderr = TRUE))
ok("--plan of cov5-mnchisq-inter-q05 names the interaction, the fit cov5-mnchisq-inter on cov5-mnchisq, and 45 minutes or the cache",
   any(grepl("Analysis: cov5-mnchisq-inter-q05 (", out, fixed = TRUE)) && any(grepl("tested ploidy by heat interaction", out, fixed = TRUE)) &&
     any(grepl("Fit of step 05: methylkit_cov5-mnchisq-inter.rds (the interaction contrast; the united table and the ploidy and heat contrasts from methylkit_cov5-mnchisq.rds, the fit of the analysis cov5-mnchisq), ", out, fixed = TRUE)) &&
     any(grepl("about 45 minutes|in the cache \\(reused", out)))
out <- suppressWarnings(system2(rscript, c("run_pipeline.R", "--plan", "--overdispersion"), stdout = TRUE, stderr = TRUE))
ok("--plan of mnchisq says contrasts of the paper and keeps its fit line",
   any(grepl("contrasts of the paper, C-to-T", out, fixed = TRUE)) && any(grepl("Fit of step 05: methylkit_mnchisq.rds, ", out, fixed = TRUE)))

# ---- Maps -------------------------------------------------------------------------------------
read_csv_file <- function(f) utils::read.csv(f, stringsAsFactors = FALSE, colClasses = "character", na.strings = NULL, check.names = FALSE,
                                              encoding = "UTF-8")
base_file <- tempfile(fileext = ".csv")
system2("git", c("show", paste0(base, ":03_input_files/manuscript/manuscript_values_map.csv")), stdout = base_file)
base_map <- read_csv_file(base_file)
map <- read_csv_file("03_input_files/manuscript/manuscript_values_map.csv")
# expression_paper_recipe (2026-11-20, after BASE) is derived from expression_legacy by
# 07_documentation/toggles_migration/paper_recipe_column.py; tests/check_paper_recipe_legacy.R checks it.
same <- vapply(setdiff(names(map), c("expression_optional", "expression_paper_recipe")), function(n) identical(map[[n]], base_map[[n]]), logical(1))
ok("values map: every column except expression_optional and expression_paper_recipe equals BASE",
   all(same) && identical(setdiff(names(base_map), names(map)), "expression_corrected") &&
     identical(setdiff(names(map), names(base_map)), c("expression_paper_recipe", "expression_optional")))
e <- map$expression_optional
e_wo_diag <- gsub('diagnostic\\(\\s*"[^"]*"\\)', "", e)
ok("values map: 263 expression_optional; no _corrected table name, no literal \"corrected\"",
   sum(nzchar(e)) == 263 && !any(grepl("_corrected", e_wo_diag)) && !any(grepl('"corrected"', e, fixed = TRUE)))
for (f in c("manuscript_values_map.csv", "manuscript_items_map.csv")) {
  raw <- readBin(file.path("03_input_files/manuscript", f), "raw", 1e7)
  n_lf <- sum(raw == as.raw(10)); n_crlf <- sum(raw[-1] == as.raw(10) & raw[-length(raw)] == as.raw(13))
  ok(sprintf("%s keeps CRLF on every line", f), n_lf > 0 && n_lf == n_crlf)
}
items <- read_csv_file("03_input_files/manuscript/manuscript_items_map.csv")
ok("items map: analyses legacy, optional, comparisons or empty",
   all(items$analysis %in% c("", "legacy", "optional", "comparisons")) && !anyDuplicated(items$destination))

res <- rbindlist(results)
cat(sprintf("\nT0: %d of %d checks pass\n", sum(res$pass), nrow(res)))
quit(status = if (all(res$pass)) 0 else 1)
