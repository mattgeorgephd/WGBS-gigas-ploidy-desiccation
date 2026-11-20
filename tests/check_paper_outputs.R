# tests/check_paper_outputs.R
#
# T1 of the settings layout (settings_build_spec.md, section 8): after a default
# run (Rscript run_pipeline.R, or a run of some of its steps, for example
# --steps 5:13 with the fit of step 05 in the cache), the paper's outputs are
# those expected. Run from the repository root:
#   Rscript tests/snapshot_paper.R snapshot <expected.csv>     before the run
#   Rscript run_pipeline.R [--steps ...]
#   Rscript tests/check_paper_outputs.R <expected.csv>         after the run
# Exit status 1 when any of the four checks fails:
#  1. Snapshot. Every file of <expected.csv> has the bytes it had before the run
#     (for .gz files the decompressed content, because the gzip header names the
#     operating system), PDFs and session_info.txt left out (their time stamps
#     change on every run), and README files (documentation, which no step
#     writes); no paper file is added or removed. The one exception
#     is 04_output/01_fetch_inputs/tables/input_inventory.csv, which lists the md5
#     of every input: its rows of the two manuscript maps, of the alignment
#     manifest and of README files may differ from the commit T1_INVENTORY_BASE
#     (default HEAD), every other row must not.
#  2. BASE (spec 7.3). Against tests/baseline/paper_manifest_BASE.csv, the
#     manifest of every tracked file at BASE (66c7121, or T1_BASE), and the files
#     of that commit: every paper file of BASE outside table 7.3 has its BASE md5
#     (PDFs, session_info.txt, the inventory and the README files, which are
#     documentation, left out); every paper file of
#     BASE that is gone was moved (a file of the same name is under a
#     not_adopted/ folder); each file of table 7.3 satisfies its projection
#     check (the BASE file minus the rows or columns of the optional analyses).
#  3. git status --porcelain of 04_output and 05_manuscript, README files
#     (documentation) left out: outside not_adopted/, no line other than PDFs,
#     session_info.txt and input_inventory.csv (whose rows check 1 tests);
#     under not_adopted/, only the lines that the snapshot recorded before the
#     run, with the same md5 (a default run writes nothing there).

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript tests/check_paper_outputs.R <expected.csv>", call. = FALSE)
fails <- character(0)
fail <- function(what) { fails <<- c(fails, what); cat("FAIL", what, "\n") }
md5_inner <- function(f) {
  tmp <- tempfile(); on.exit(unlink(tmp))
  con <- gzfile(f, "rb"); out <- file(tmp, "wb")
  repeat { b <- readBin(con, "raw", 1e7); if (length(b) == 0) break; writeBin(b, out) }
  close(con); close(out)
  unname(tools::md5sum(tmp))
}
md5_now <- function(f) if (endsWith(f, ".gz")) md5_inner(f) else unname(tools::md5sum(f))
git_show <- function(commit, path) {
  f <- tempfile(fileext = paste0(".", tools::file_ext(path)))
  status <- system2("git", c("show", paste0(commit, ":", path)), stdout = f, stderr = FALSE)
  if (!identical(as.integer(status), 0L)) stop("git show ", commit, ":", path, " failed", call. = FALSE)
  f
}
skip_file <- function(p) grepl("\\.pdf$", p) | basename(p) == "session_info.txt"
inventory <- "04_output/01_fetch_inputs/tables/input_inventory.csv"
now <- c(list.files("04_output", recursive = TRUE, full.names = TRUE), list.files("05_manuscript", recursive = TRUE, full.names = TRUE))
now <- now[!grepl("(^|/)not_adopted(/|$)", now) & !grepl("/cache/", now)]

# ---- 1. Snapshot ---------------------------------------------------------------------------------
expected <- fread(args[1], colClasses = "character")[md5 != "folder"]
paper <- expected$path[!startsWith(expected$path, "06_diagnostics/")]
added <- setdiff(now, paper); removed <- setdiff(paper, now)
skip <- skip_file(expected$path)
check <- expected[!skip & path != inventory & basename(path) != "README.md" & file.exists(path)]
check[, md5_after := vapply(path, md5_now, "")]
check[, same := ifelse(endsWith(path, ".gz"), md5_content == md5_after, md5 == md5_after)]
cat(sprintf("1. snapshot: %d paper files compared byte for byte (PDFs and session_info.txt left out: %d; README files too); %d differ; %d added, %d removed\n",
            nrow(check), sum(skip), sum(!check$same), length(added), length(removed)))
if (any(!check$same)) { fail("files differ from the snapshot"); print(check[same == FALSE, .(path)]) }
if (length(added)) { fail("paper files added"); cat("  ", added, sep = "\n  ") }
if (length(removed)) { fail("paper files removed"); cat("  ", removed, sep = "\n  ") }
inventory_ok <- TRUE
if (inventory %in% expected$path) {
  inv_base <- Sys.getenv("T1_INVENTORY_BASE", "HEAD")
  old <- fread(git_show(inv_base, inventory), colClasses = "character"); new <- fread(inventory, colClasses = "character")
  key <- c("location", "file")
  row_text <- function(d) setNames(do.call(paste, c(d[, setdiff(names(d), key), with = FALSE], sep = "\t")), paste(d$location, d$file))
  o <- row_text(old); n <- row_text(new)
  rows <- union(names(o), names(n))
  differs <- rows[is.na(o[rows]) | is.na(n[rows]) | o[rows] != n[rows]]
  allowed <- grepl("manuscript_(values|items)_map\\.csv$|gannet_manifest_alignments\\.csv$|README\\.md$", differs)
  inventory_ok <- all(allowed)
  cat(sprintf("   input_inventory.csv: %d rows; %d differ from %s: %s\n", nrow(new), length(differs), inv_base, paste(differs, collapse = ", ")))
  if (!inventory_ok) fail(paste("input_inventory.csv rows that may not change:", paste(differs[!allowed], collapse = ", ")))
}

# ---- 2. BASE and the projections of spec 7.3 --------------------------------------------------------
base_commit <- Sys.getenv("T1_BASE", "66c7121")
base <- fread("tests/baseline/paper_manifest_BASE.csv", colClasses = "character")
base <- base[(startsWith(path, "04_output/") | startsWith(path, "05_manuscript/")) & !grepl("(^|/)not_adopted(/|$)", path)]
lines_of <- function(f) readLines(f, warn = FALSE, encoding = "UTF-8")
text_of <- function(f) fread(f, colClasses = "character", na.strings = NULL)
# 7.3, rows of the optional analyses removed: header and every line equal, in order, to the BASE lines with mode == "legacy".
legacy_rows <- function(p) {
  b_file <- git_show(base_commit, p); b <- lines_of(b_file); bt <- text_of(b_file)
  if (nrow(bt) != length(b) - 1) stop("a field of ", p, " spans lines")
  identical(lines_of(p), c(b[1], b[-1][bt$mode == "legacy"]))
}
# 7.3, columns of the optional analyses removed: exactly `cols` gone, every other field equal to BASE as text.
columns_removed <- function(p, cols) {
  b <- text_of(git_show(base_commit, p)); n <- text_of(p)
  setequal(setdiff(names(b), names(n)), cols) && length(setdiff(names(n), names(b))) == 0 &&
    isTRUE(all.equal(b[, names(n), with = FALSE], n, check.attributes = FALSE))
}
# 7.3, rows of the gene-level part removed: same header, the remaining lines are BASE lines in their order.
subsequence <- function(p) {
  b <- lines_of(git_show(base_commit, p)); n <- lines_of(p)
  identical(b[1], n[1]) && !anyNA(m <- match(n[-1], b[-1])) && !is.unsorted(m, strictly = TRUE)
}
items_rows <- function(p) {
  b <- text_of(git_show(base_commit, p)); n <- text_of(p)
  keep <- b[analysis %in% c("", "legacy") & !grepl("10_go_enrichment_gene_level", source)]
  identical(names(b), names(n)) && isTRUE(all.equal(keep, n, check.attributes = FALSE))
}
projected <- list(
  "04_output/04_cpg_coverage/tables/cpg_counts_by_coverage.csv"              = function(p) legacy_rows(p),
  "04_output/04_cpg_coverage/tables/cpg_counts_by_feature_track.csv"         = function(p) legacy_rows(p),
  "04_output/04_cpg_coverage/tables/methylation_distribution_by_library.csv" = function(p) legacy_rows(p),
  "04_output/04_cpg_coverage/tables/methylation_distribution_by_group.csv"   = function(p) legacy_rows(p),
  "04_output/04_cpg_coverage/tables/library_correlation_summary.csv"         = function(p) legacy_rows(p),
  "04_output/04_cpg_coverage/tables/pca_variance.csv"                        = function(p) legacy_rows(p),
  "04_output/05_differential_methylation/tables/dml_counts_by_mode.csv"      = function(p) legacy_rows(p),
  "04_output/08_gene_annotation/tables/annotation_summary.csv"               = function(p) legacy_rows(p),
  "04_output/11_semantic_similarity/tables/named_terms_results_3_4.csv"      = function(p) columns_removed(p, c(
    "gene_level_genes_annotated", "gene_level_genes_with_dml", "gene_level_p_wallenius", "gene_level_q_wallenius", "gene_level_result")),
  "04_output/11_semantic_similarity/tables/semantic_summary.csv"             = function(p) subsequence(p),
  "04_output/11_semantic_similarity/tables/semantic_terms_dropped.csv"       = function(p) subsequence(p),
  "05_manuscript/manuscript_values.csv"         = function(p) columns_removed(p, c("regenerated_corrected", "corrected_differs", "expression_corrected")),
  "05_manuscript/manuscript_values_summary.csv" = function(p) columns_removed(p, "change_in_corrected_mode"),
  "04_output/13_manuscript_values/tables/expression_errors.csv"        = function(p) columns_removed(p, "expression_corrected"),
  "04_output/13_manuscript_values/tables/manuscript_items_copied.csv" = function(p) items_rows(p))
gone <- base$path[!file.exists(base$path)]
moved_names <- basename(c(list.files("04_output/not_adopted", recursive = TRUE), list.files("05_manuscript/not_adopted", recursive = TRUE)))
not_moved <- gone[!basename(gone) %in% moved_names]
same_base <- base[file.exists(path) & !skip_file(path) & path != inventory & basename(path) != "README.md" & !path %in% names(projected)]
same_base[, md5_after := vapply(path, md5_now, "")]
same_base[, same := ifelse(endsWith(path, ".gz"), md5_content == md5_after, md5 == md5_after)]
proj <- vapply(names(projected), function(p) file.exists(p) && isTRUE(tryCatch(projected[[p]](p), error = function(e) FALSE)), logical(1))
cat(sprintf("2. BASE %s: %d paper files of BASE compared with their BASE md5, %d differ; %d gone, of which %d not under a not_adopted/ folder; %d of %d projections of spec 7.3 hold\n",
            base_commit, nrow(same_base), sum(!same_base$same), length(gone), length(not_moved), sum(proj), length(proj)))
if (any(!same_base$same)) { fail("paper files differ from BASE"); print(same_base[same == FALSE, .(path)]) }
if (length(not_moved)) { fail("paper files of BASE gone and not moved"); cat("  ", not_moved, sep = "\n  ") }
if (!all(proj)) { fail("projections of spec 7.3 that do not hold"); cat("  ", names(proj)[!proj], sep = "\n  ") }

# ---- 3. git status ---------------------------------------------------------------------------------------
porcelain <- function() system2("git", c("status", "--porcelain", "--untracked-files=all", "--", "04_output", "05_manuscript"), stdout = TRUE)
st <- porcelain()
st_path <- sub(".* -> ", "", substring(st, 4))
under_na <- grepl("(^|/)not_adopted/", st_path)
docs <- basename(st_path) == "README.md"          # documentation, not an output of a step
outside <- st[!under_na & !docs & !skip_file(st_path) & !(st_path == inventory & inventory_ok)]
cat(sprintf("3. git status of 04_output and 05_manuscript: %d lines; outside not_adopted/ (PDFs, session_info.txt, README files and a checked inventory left out): %d\n",
            length(st), length(outside)))
if (length(outside)) { fail("git status lines outside not_adopted/"); cat("  ", outside, sep = "\n  ") }
before_file <- paste0(args[1], ".not_adopted_status.csv")
before <- if (file.exists(before_file)) fread(before_file, colClasses = "character", strip.white = FALSE) else data.table(line = character(0), md5 = character(0))
na_now <- data.table(line = st[under_na & !docs], path = st_path[under_na & !docs])
na_now[, md5 := ifelse(file.exists(path), unname(tools::md5sum(path)), "absent")]
unexpected <- na_now[!paste(line, md5) %in% paste(before$line, before$md5), line]
cat(sprintf("   under not_adopted/: %d lines, %d recorded before the run (%s); %d not expected\n", nrow(na_now), nrow(before),
            if (file.exists(before_file)) basename(before_file) else "no record: none expected", length(unexpected)))
if (length(unexpected)) { fail("git status lines under not_adopted/ that the run made"); cat("  ", unexpected, sep = "\n  ") }

cat(if (length(fails)) sprintf("\nT1: %d check(s) failed\n", length(fails)) else "\nT1: all checks pass\n")
quit(status = if (length(fails)) 1 else 0)
