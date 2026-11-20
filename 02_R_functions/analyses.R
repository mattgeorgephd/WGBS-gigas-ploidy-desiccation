# analyses.R
# The analyses of a run: the paper's and the optional ones (NOT ADOPTED).
#
# One run computes ONE analysis. Its id is made of the tokens of the optional
# settings that are on, in this fixed order, joined by "-":
#   corrected  the coverage filter of the Methods (99.9th percentile; every line read)
#   cov5       5 reads per CpG in each library instead of 10
#   mnchisq    methylKit's overdispersion correction, overdispersion = "MN", test = "Chisq"
#   ctexcl     CpGs at C-to-T variants removed (Venkataraman et al. 2022)
# With no token on, the id is "legacy": the analysis of the paper. The paper's
# outputs are in 04_output/<step>/ and 05_manuscript/; the outputs of any other
# analysis are in 04_output/not_adopted/<id>/ and 05_manuscript/not_adopted/.
# The settings and switches are described in run_config.R and run_pipeline.R;
# the layout in 04_output/not_adopted/README.md.
#
# run_config.R sources this file in its checks, before the rest of
# 02_R_functions/ is sourced, so that analysis_spec() can check WGBS_ANALYSIS.
# Every other function here is used by the steps, after config exists.

analysis_tokens <- c("corrected", "cov5", "mnchisq", "ctexcl")

# The id of a set of tokens, in the fixed order; "legacy" when there is none.
analysis_id <- function(tokens) {
  t <- analysis_tokens[analysis_tokens %in% tokens]
  if (length(t) == 0) "legacy" else paste(t, collapse = "-")
}

# What an analysis id means. Stops on an id that is not one of the 16 valid ids.
analysis_spec <- function(id) {
  valid <- "an analysis id: \"legacy\" (the paper) or the tokens corrected, cov5, mnchisq, ctexcl (in this order, each at most once) joined by \"-\""
  if (!is.character(id) || length(id) != 1 || is.na(id) || !nzchar(id)) stop("Not ", valid, ".", call. = FALSE)
  tokens <- if (id == "legacy") character(0) else strsplit(id, "-", fixed = TRUE)[[1]]
  if (!all(tokens %in% analysis_tokens) || anyDuplicated(tokens) || !identical(analysis_id(tokens), id)) {
    stop("'", id, "' is not ", valid, "; for example cov5, corrected-cov5, cov5-mnchisq-ctexcl.", call. = FALSE)
  }
  mn <- "mnchisq" %in% tokens
  list(id = id, paper = id == "legacy", tokens = tokens,
       filter = if ("corrected" %in% tokens) "corrected" else "legacy",
       min_coverage = if ("cov5" %in% tokens) 5 else 10,
       overdispersion = if (mn) "MN" else "none",
       test = if (mn) "Chisq" else "not passed",
       ct_excluded = "ctexcl" %in% tokens,
       fit_id = analysis_id(setdiff(tokens, "ctexcl")))
}

# Analyses in the order of every table that lists several: the paper first, then
# the other ids in byte order.
analysis_order <- function(ids) {
  ids <- unique(ids)
  c(intersect("legacy", ids), sort(setdiff(ids, "legacy"), method = "radix"))
}

# ---- Folders -----------------------------------------------------------------------
not_adopted_path <- function(...) here::here("04_output", "not_adopted", ...)

# Where an analysis writes the outputs of a step. The paper's analysis in a paper
# run: the paper folder 04_output/<step>/<sub>/. Any analysis in an optional run
# (the paper's included, when only the gene-level GO or the diagnostics are on):
# 04_output/not_adopted/<id>/<step>/<sub>/. Created if absent.
analysis_dir <- function(a, step, sub = NULL) {
  if (identical(config$run_kind, "paper")) {
    if (!a$paper) stop("A paper run writes the paper's analysis only, not '", a$id, "'.", call. = FALSE)
    return(output_dir(step, sub))
  }
  d <- if (is.null(sub)) not_adopted_path(a$id, step) else not_adopted_path(a$id, step, sub)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
analysis_path <- function(a, step, sub, file) file.path(analysis_dir(a, step, sub), file)

# Where the outputs of a step of an analysis are read from (never creates a
# folder). The paper's analysis: the paper folder, or not_adopted/legacy/ for a
# file that only an optional run writes for it (steps 10 and 11, gene level).
analysis_input <- function(a, step, sub, file) {
  id <- if (is.list(a)) a$id else a
  if (id == "legacy") {
    p <- paper_path(step, sub, file)
    if (file.exists(p)) return(p)
  }
  not_adopted_path(id, step, sub, file)
}

# ---- Stamp of an analysis folder ---------------------------------------------------
# 04_output/not_adopted/<id>/analysis.csv, one row: the settings of the analysis,
# the md5 of the fingerprint of its fit (step 05), the md5 of the variant list
# (ctexcl), the commit, the time and the steps done. Step 04 creates it (and the
# first step of the paper's analysis in an optional run); each later step adds
# itself to steps_done. A later step stops if the stamp's settings differ from
# the analysis of the run or a step it reads from is not done.
stamp_columns <- c("id", "filter", "min_coverage", "overdispersion", "test", "ct_excluded", "fit_id")
stamp_file <- function(id) not_adopted_path(id, "analysis.csv")
read_stamp <- function(id) {
  f <- stamp_file(id)
  if (!file.exists(f)) return(NULL)
  data.table::fread(f, colClasses = "character")
}
git_commit <- function() {
  out <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(here::here()), "rev-parse", "--short", "HEAD"),
                                           stdout = TRUE, stderr = FALSE)), error = function(e) character(0))
  if (length(out) == 1 && nzchar(out)) out else "unknown"
}
write_stamp <- function(a, step, fit_fingerprint_md5 = NULL, ct_variant_list_md5 = NULL) {
  if (identical(config$run_kind, "paper")) return(invisible(NULL))     # the paper folder has no stamp
  number <- substr(step, 1, 2)
  old <- read_stamp(a$id)
  new <- old
  if (is.null(old) || number == "04") {
    new <- data.table::data.table(id = a$id, filter = a$filter, min_coverage = format(a$min_coverage),
                                  overdispersion = a$overdispersion, test = a$test, ct_excluded = as.character(a$ct_excluded),
                                  fit_id = a$fit_id, fit_fingerprint_md5 = "", ct_variant_list_md5 = "",
                                  git_commit = "", written_utc = "", steps_done = "")
  }
  done <- setdiff(strsplit(new$steps_done, ";", fixed = TRUE)[[1]], "")
  new$steps_done <- paste(sort(unique(c(done, number))), collapse = ";")
  if (!is.null(fit_fingerprint_md5)) new$fit_fingerprint_md5 <- fit_fingerprint_md5
  if (!is.null(ct_variant_list_md5)) new$ct_variant_list_md5 <- ct_variant_list_md5
  new$git_commit <- git_commit()
  new$written_utc <- format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC")
  dir.create(not_adopted_path(a$id), recursive = TRUE, showWarnings = FALSE)
  write_table(new, stamp_file(a$id))
  invisible(new)
}
# Stops unless the folder of `a` holds a stamp of the same analysis that lists
# every step of `needs` (step folder names or two-digit numbers) as done. With
# per_cpg = TRUE (steps 06 to 08, which read united_cpg_<id>.csv.gz) it also
# stops when the united table of step 05 is absent: the per-CpG tables of an
# optional analysis other than 'corrected' are not tracked by git (.gitignore),
# so a clone holds the stamp and the small tables but not those. Not used for
# the paper's analysis: its inputs are the paper's outputs.
check_stamp <- function(a, needs = character(0), per_cpg = FALSE) {
  if (a$paper) return(invisible(TRUE))
  s <- read_stamp(a$id)
  if (is.null(s)) {
    stop("04_output/not_adopted/", a$id, "/analysis.csv is absent: step 04 has not run for the analysis '", a$id,
         "'. Run Rscript run_pipeline.R with its switches (or WGBS_ANALYSIS=", a$id, " for the steps from 04).", call. = FALSE)
  }
  want <- c(id = a$id, filter = a$filter, min_coverage = format(a$min_coverage), overdispersion = a$overdispersion,
            test = a$test, ct_excluded = as.character(a$ct_excluded), fit_id = a$fit_id)
  differs <- names(want)[unlist(s[, names(want), with = FALSE]) != want]
  if (length(differs) > 0) {
    stop("The stamp 04_output/not_adopted/", a$id, "/analysis.csv does not describe the analysis '", a$id, "' (", paste(differs, collapse = ", "),
         " differ). Run step 04 for it again.", call. = FALSE)
  }
  done <- strsplit(s$steps_done, ";", fixed = TRUE)[[1]]
  missing <- setdiff(substr(needs, 1, 2), done)
  if (length(missing) > 0) {
    stop("Step(s) ", paste(missing, collapse = ", "), " have not run for the analysis '", a$id, "' (stamp: steps ",
         s$steps_done, "). Run them first.", call. = FALSE)
  }
  if (per_cpg) require_united_table(a)
  invisible(TRUE)
}
# The table of every tested CpG of an analysis (step 05), which steps 06 to 08
# read. Stops with the command that writes it again when it is absent.
united_table_file <- function(a) not_adopted_path(a$id, "05_differential_methylation", "tables", sprintf("united_cpg_%s.csv.gz", a$id))
require_united_table <- function(a) {
  if (a$paper || file.exists(united_table_file(a))) return(invisible(TRUE))
  stop("04_output/not_adopted/", a$id, "/05_differential_methylation/tables/united_cpg_", a$id, ".csv.gz is absent: ",
       "rerun step 05 for ", a$id, " (per-CpG tables are not tracked): WGBS_ANALYSIS=", a$id,
       " Rscript 01_pipeline/05_differential_methylation.R (half a minute with the fit in the cache, else a new fit).", call. = FALSE)
}
# Whether the stamp of the analysis folder `id` lists `step` (a step folder name
# or a two-digit number) as done.
stamp_lists_step <- function(id, step) {
  s <- read_stamp(id)
  !is.null(s) && substr(step, 1, 2) %in% strsplit(s$steps_done, ";", fixed = TRUE)[[1]]
}
# Step 04 of an optional analysis starts its folder afresh: it removes the
# stamp and the folders of steps 05 and later of the analysis (and the copies
# of its items under 05_manuscript/not_adopted/<id>/) before it writes, so that
# no table of an earlier run of a later step survives a rerun that stops
# halfway. The fits of step 05 are in the cache folder of the paper's step 05
# and are kept. Returns the folders removed.
clear_analysis_folder <- function(a) {
  if (!identical(config$run_kind, "optional") || a$paper) return(invisible(character(0)))
  d <- list.dirs(not_adopted_path(a$id), recursive = FALSE, full.names = FALSE)
  later <- d[grepl("^[0-9]{2}_", d) & as.integer(substr(d, 1, 2)) >= 5]
  unlink(not_adopted_path(a$id, later), recursive = TRUE)
  unlink(stamp_file(a$id))
  unlink(manuscript_path("not_adopted", a$id), recursive = TRUE)
  invisible(later)
}
# Ids of the analysis folders under 04_output/not_adopted/ that hold a stamp.
analysis_folders <- function() {
  d <- list.dirs(not_adopted_path(), recursive = FALSE, full.names = FALSE)
  d <- d[file.exists(not_adopted_path(d, "analysis.csv"))]
  analysis_order(d)
}

# ---- Comparison tables (04_output/not_adopted/comparisons/) -------------------------
# The table `file` of `step` of every analysis on disk, side by side, the paper
# first: the paper's copy (paper folder, or not_adopted/legacy/), then each
# analysis folder in analysis order whose stamp lists `step` as done. Lines are copied as they are, so every row
# keeps its bytes; `arrange` (a function of the bound table, read with fread)
# returns the order of the rows where the paper's table interleaves analyses.
# Writes comparisons/<step>/tables/<file> and comparisons/analyses_included.csv.
rebuild_comparison <- function(step, file, arrange = NULL) {
  if (!identical(config$run_kind, "optional")) stop("rebuild_comparison() runs in optional runs only.", call. = FALSE)
  ids <- analysis_order(c("legacy", analysis_folders()))
  paths <- vapply(ids, function(id) analysis_input(id, step, "tables", file), "")
  # A table in an analysis folder counts only while the stamp lists its step:
  # a table left by a run that stopped, or by an earlier run, is not copied.
  usable <- file.exists(paths) & (paths == paper_path(step, "tables", file) | vapply(ids, stamp_lists_step, logical(1), step = step))
  ids <- ids[usable]; paths <- paths[usable]
  if (length(paths) == 0) return(invisible(NULL))
  header <- NULL; body <- character(0); tables <- list()
  for (i in seq_along(paths)) {
    l <- readLines(paths[i], warn = FALSE)
    if (is.null(header)) header <- l[1]
    if (!identical(l[1], header)) stop("The columns of ", paths[i], " differ from those of ", paths[1], ".", call. = FALSE)
    body <- c(body, l[-1])
    if (!is.null(arrange)) tables[[i]] <- data.table::fread(paths[i], colClasses = "character")
  }
  if (!is.null(arrange) && length(body) > 0) {
    bound <- data.table::rbindlist(tables)
    if (nrow(bound) != length(body)) stop("Cannot order ", file, ": a field spans lines.", call. = FALSE)
    body <- body[arrange(bound)]
  }
  out <- not_adopted_path("comparisons", step, "tables")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  con <- file(file.path(out, file), "wb"); writeLines(c(header, body), con, sep = "\n"); close(con)
  write_analyses_included()
  invisible(file.path(out, file))
}
write_analyses_included <- function() {
  ids <- analysis_folders()
  rows <- lapply(ids, read_stamp)
  paper <- data.table::data.table(id = "legacy", filter = "legacy", min_coverage = "10", overdispersion = "none", test = "not passed",
                                  ct_excluded = "FALSE", fit_id = "legacy", fit_fingerprint_md5 = "", ct_variant_list_md5 = "",
                                  git_commit = "", written_utc = "", steps_done = "the paper's outputs (04_output/<step>/)")
  if ("legacy" %in% ids) {
    paper <- rows[[match("legacy", ids)]]
    paper$steps_done <- paste0("the paper's outputs (04_output/<step>/); in not_adopted/legacy/: steps ", paper$steps_done)
  }
  tab <- data.table::rbindlist(c(list(paper), rows[!ids %in% "legacy"]), fill = TRUE)
  dir.create(not_adopted_path("comparisons"), recursive = TRUE, showWarnings = FALSE)
  write_table(tab, not_adopted_path("comparisons", "analyses_included.csv"))
}

# ---- C-to-T variant list (option ctexcl) --------------------------------------------
ct_variant_files <- function() {
  d <- not_adopted_path("ct_variants", "tables")
  c(cpgs = file.path(d, "ct_variant_cpgs.csv.gz"), assessability = file.path(d, "ct_variant_assessability.csv.gz"),
    index = file.path(d, "ct_variant_cpg_index.csv"))
}
ct_variant_build_command <- paste(
  "WGBS_CACHE_DIR=<download cache> Rscript 01_pipeline/optional/ct_variant_loci.R   (from the repository root; downloads the",
  "22 alignments one at a time, 109 GB in all, about 7 to 9 hours; see 01_pipeline/optional/README.md)")
ct_variant_list_present <- function() all(file.exists(ct_variant_files()))
# The CpGs listed by the 2022 rule (column `listed`): chr, start (0-based C, as
# in the coverage files). Stops with the command that makes the list when the
# variant step has not finished.
ct_variant_table <- function() {
  f <- ct_variant_files()
  if (!all(file.exists(f))) {
    stop("The option --exclude-ct-variants needs the C-to-T variant list, and ", paste(basename(f[!file.exists(f)]), collapse = ", "),
         " is absent from 04_output/not_adopted/ct_variants/tables/. Make it first with:\n  ", ct_variant_build_command, call. = FALSE)
  }
  v <- data.table::fread(f[["cpgs"]], select = c("chr", "start", "listed"),
                         colClasses = list(character = "chr", integer = "start", character = "listed"))
  v <- v[toupper(v$listed) == "TRUE", c("chr", "start")]
  data.table::setkeyv(v, c("chr", "start"))
  v
}
ct_variant_list_md5 <- function() unname(tools::md5sum(ct_variant_files()[["cpgs"]]))
# Every CpG of the reference with the evidence of the variant screen:
# crick_pooled (Crick C + T reads at the C over the 22 libraries) and
# assessable (crick_pooled at least config$ct_variants$min_crick_assessable,
# the fewest evidence reads with which any call set can list the site). The rows
# of ct_variant_assessability.csv.gz carry no coordinates; they are the plus-strand
# "CG" of the reference, sequences in FASTA order, positions increasing, rebuilt
# here and checked against ct_variant_cpg_index.csv per sequence.
ct_variant_assessability <- function() {
  f <- ct_variant_files()
  ct_variant_table()                                   # stops with the command when the files are absent
  index <- data.table::fread(f[["index"]])
  genome <- Biostrings::readDNAStringSet(fetch_gannet_file("cgigas_uk_roslin_v1_genomic-mito.fa", config))
  names(genome) <- sub(" .*", "", names(genome))
  coords <- data.table::rbindlist(lapply(seq_along(genome), function(j) {
    data.table::data.table(chr = names(genome)[j], start = BiocGenerics::start(Biostrings::matchPattern("CG", genome[[j]])) - 1L)
  }))
  md5_text <- function(starts) { tmp <- tempfile(); on.exit(unlink(tmp)); writeLines(format(starts, scientific = FALSE, trim = TRUE), tmp); unname(tools::md5sum(tmp)) }
  no_cg <- setdiff(names(genome), unique(coords$chr))
  rm(genome); invisible(gc(verbose = FALSE))
  # A sequence without a "CG" has no row in the assessability table; the index
  # of the variant step would list it with n_cpg 0 and no position. The reference
  # of the study has none (237 sequences, each with CpGs); stop clearly if one appears.
  if (length(no_cg) > 0 || any(index$n_cpg == 0)) {
    stop("The reference holds sequence(s) without a CpG (", paste(head(union(no_cg, index$chr[index$n_cpg == 0]), 5), collapse = ", "),
         "); ct_variant_assessability() places rows only on sequences with CpGs and cannot read this index.", call. = FALSE)
  }
  per_chr <- coords[, .(n_cpg = .N, first_start = start[1], last_start = start[.N], md5_starts = md5_text(start)), by = chr]
  per_chr[, first_row := cumsum(c(1L, head(n_cpg, -1L)))]          # the rows are the sequences in FASTA order
  check <- merge(index[, .(chr, n_cpg, first_row, first_start, last_start, md5_starts)], per_chr, by = "chr", suffixes = c("_index", "_genome"))
  if (nrow(check) != nrow(index) || nrow(index) != nrow(per_chr) ||
      check[, any(n_cpg_index != n_cpg_genome | md5_starts_index != md5_starts_genome)]) {
    stop("The rows of ct_variant_assessability.csv.gz cannot be placed: the CpGs of the reference differ from ct_variant_cpg_index.csv.", call. = FALSE)
  }
  if (!identical(index$chr, per_chr$chr) || check[, any(first_row_index != first_row_genome)]) {
    stop("The rows of ct_variant_assessability.csv.gz cannot be placed: the order of the sequences in ct_variant_cpg_index.csv ",
         "(column first_row) is not the order of the reference.", call. = FALSE)
  }
  a <- data.table::fread(f[["assessability"]], select = "crick_pooled")
  if (nrow(a) != nrow(coords)) stop("ct_variant_assessability.csv.gz has ", nrow(a), " rows; the reference has ", nrow(coords), " CpGs.", call. = FALSE)
  coords[, crick_pooled := a$crick_pooled]
  coords[, assessable := crick_pooled >= config$ct_variants$min_crick_assessable]
  data.table::setkeyv(coords, c("chr", "start"))
  coords
}

# ---- Run kinds -----------------------------------------------------------------------
# End a step that has nothing to do in this run kind. `kind` is the run kind the
# step needs ("paper" or "optional"); `step_does` says what the step does, and
# opens the message. Under Rscript the step ends here with exit status 0, before
# it reads or writes anything, so that a job script that runs steps one by one
# goes on; a script sourced in an interactive session stops with an error
# (quit() would close the session). With kind "optional_gene_level" the step
# needs an optional run with the gene-level GO enrichment on.
require_run_kind <- function(kind, step_does) {
  ok <- switch(kind,
               paper = identical(config$run_kind, "paper"),
               optional = identical(config$run_kind, "optional"),
               optional_gene_level = identical(config$run_kind, "optional") && isTRUE(config$addons$gene_level_go),
               stop("Unknown run kind: ", kind))
  if (ok) return(invisible(TRUE))
  wanted <- switch(kind, paper = "the paper run (no optional setting on)", optional = "an optional run",
                   optional_gene_level = "an optional run with the gene-level GO enrichment on (--gene-level-go)")
  not_run <- sprintf("%s, and is part of %s; this is %s: nothing is run and no file is written",
                     step_does, wanted, run_kind_text())
  if (interactive()) stop(not_run, call. = FALSE)
  log_msg("%s", not_run)
  quit(save = "no", status = 0)
}
run_kind_text <- function() {
  if (identical(config$run_kind, "paper")) return("the paper run")
  addons <- c(if (isTRUE(config$addons$gene_level_go)) "gene-level GO",
              if (length(config$addons$diagnostics) > 0) paste0("diagnostics (", paste(config$addons$diagnostics, collapse = ","), ")"))
  sprintf("an optional run (analysis '%s'%s)", config$analysis$id,
          if (length(addons)) paste0("; ", paste(addons, collapse = "; ")) else "")
}

# ---- Diagnostics -----------------------------------------------------------------------
# The diagnostics of 06_diagnostics/ compare the paper's analysis ("legacy") with
# the coverage filter of the Methods ("corrected"), whose outputs are in
# 04_output/not_adopted/corrected/. Both must be present.
diagnostic_modes <- function() {
  modes <- c("legacy", "corrected")
  if (!file.exists(not_adopted_path("corrected", "analysis.csv"))) {
    stop("The diagnostics read the outputs of the analysis 'corrected', and 04_output/not_adopted/corrected/ holds none. ",
         "Make them with: Rscript run_pipeline.R --coverage-filter", call. = FALSE)
  }
  modes
}
# Path of a pipeline output read by a diagnostic (never creates a folder): a file
# of the analysis 'corrected' (name ending in _corrected) from
# 04_output/not_adopted/corrected/, every other file from the paper folder.
diagnostic_input <- function(step, sub, file) {
  if (grepl("_corrected\\.[A-Za-z.]+$", file)) not_adopted_path("corrected", step, sub, file) else paper_path(step, sub, file)
}
