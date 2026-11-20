# tests/check_paper_outputs.R
#
# T1 of the settings layout (settings_build_spec.md, section 8): after a default
# run (Rscript run_pipeline.R), every paper output has the bytes it had before.
#   Rscript tests/snapshot_paper.R snapshot <expected.csv>     before the run
#   Rscript run_pipeline.R
#   Rscript tests/check_paper_outputs.R <expected.csv>         after the run
# Compared: md5 of every file of the snapshot (for .gz files the md5 of the
# decompressed content, because the gzip header names the operating system),
# except PDFs and session_info.txt, whose time stamps change on every run.
# 04_output/01_fetch_inputs/tables/input_inventory.csv lists the md5 of every
# input; its rows of the two manuscript maps (rewritten by the settings layout)
# and of README files may change, every other row must not. It is compared with
# the commit T1_INVENTORY_BASE (default HEAD). Also: no paper file added or removed,
# and no file of the paper run under a not_adopted/ folder changed by the run
# (git status). Exit status 1 when anything differs.

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
expected <- fread(args[1], colClasses = "character")[md5 != "folder"]
md5_inner <- function(f) {
  tmp <- tempfile(); on.exit(unlink(tmp))
  con <- gzfile(f, "rb"); out <- file(tmp, "wb")
  repeat { b <- readBin(con, "raw", 1e7); if (length(b) == 0) break; writeBin(b, out) }
  close(con); close(out)
  unname(tools::md5sum(tmp))
}
now <- c(list.files("04_output", recursive = TRUE, full.names = TRUE), list.files("05_manuscript", recursive = TRUE, full.names = TRUE))
now <- now[!grepl("(^|/)not_adopted(/|$)", now) & !grepl("/cache/", now)]
paper <- expected$path[!startsWith(expected$path, "06_diagnostics/")]
added <- setdiff(now, paper); removed <- setdiff(paper, now)
skip <- grepl("\\.pdf$", expected$path) | basename(expected$path) == "session_info.txt"
inventory <- "04_output/01_fetch_inputs/tables/input_inventory.csv"
check <- expected[!skip & path != inventory & file.exists(path)]
check[, md5_now := ifelse(endsWith(path, ".gz"), vapply(path, md5_inner, ""), unname(tools::md5sum(path)))]
check[, same := ifelse(endsWith(path, ".gz"), md5_content == md5_now, md5 == md5_now)]
cat(sprintf("T1: %d paper files compared byte for byte (PDFs and session_info.txt left out: %d); %d differ; %d added, %d removed\n",
            nrow(check), sum(skip), sum(!check$same), length(added), length(removed)))
if (any(!check$same)) print(check[same == FALSE, .(path)])
if (length(added)) cat("added:", added, sep = "\n  ")
if (length(removed)) cat("removed:", removed, sep = "\n  ")

# input_inventory.csv: every row except those of the two maps and of README files edited with the layout
inv_ok <- TRUE
if (inventory %in% expected$path) {
  old_file <- tempfile(fileext = ".csv")
  system2("git", c("show", paste0(Sys.getenv("T1_INVENTORY_BASE", "HEAD"), ":", inventory)), stdout = old_file)
  old <- fread(old_file, colClasses = "character"); new <- fread(inventory, colClasses = "character")
  key <- c("location", "file")
  row_text <- function(d) setNames(do.call(paste, c(d[, setdiff(names(d), key), with = FALSE], sep = "\t")), paste(d$location, d$file))
  o <- row_text(old); n <- row_text(new)
  paths <- union(names(o), names(n))
  differs <- paths[is.na(o[paths]) | is.na(n[paths]) | o[paths] != n[paths]]
  allowed <- grepl("manuscript_(values|items)_map\\.csv$|README\\.md$", differs)
  inv_ok <- all(allowed)
  cat(sprintf("input_inventory.csv: %d rows; %d rows differ from %s: %s\n", nrow(new), length(differs),
              Sys.getenv("T1_INVENTORY_BASE", "HEAD"), paste(differs, collapse = ", ")))
  if (!inv_ok) cat("rows that may not change:", differs[!allowed], "\n")
}
st <- system2("git", c("status", "--porcelain", "--", "04_output", "05_manuscript"), stdout = TRUE)
st <- st[!grepl("session_info\\.txt$|\\.pdf$", st)]
cat(sprintf("git status of 04_output and 05_manuscript (PDFs and session_info.txt left out): %d lines\n", length(st)))
if (length(st)) cat(st, sep = "\n")
quit(status = if (all(check$same) && !length(added) && !length(removed) && inv_ok) 0 else 1)
