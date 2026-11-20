# tests/snapshot_paper.R
#
# Manifests of output files, for the tests of the settings layout
# (07_documentation/agent_workflows/scope_inputs/settings_build_spec.md, section 8).
#
#   Rscript tests/snapshot_paper.R base <out.csv>
#       every tracked file under 04_output/, 05_manuscript/ and 06_diagnostics/:
#       path, bytes, md5, and for .gz files the md5 of the decompressed content
#       (the gzip header names the operating system). Written once, at BASE.
#   Rscript tests/snapshot_paper.R snapshot <out.csv>
#       T2: every file of the paper (under 04_output/ and 05_manuscript/ outside
#       not_adopted/ and outside the cache folders), plus the two diagnostic tables
#       that paper values read: path, bytes, md5, modification time.
#   Rscript tests/snapshot_paper.R compare <before.csv> <after.csv>
#       T2: lists every file added, removed or changed (md5 or modification time);
#       exit status 1 when anything differs.
#
# Run from the repository root. R only; reads files, writes only <out.csv>.

args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(data.table))
md5_of <- function(f) unname(tools::md5sum(f))
md5_inner <- function(f) {
  if (!grepl("\\.gz$", f)) return(NA_character_)
  tmp <- tempfile(); on.exit(unlink(tmp))
  con <- gzfile(f, "rb"); out <- file(tmp, "wb")
  repeat { b <- readBin(con, "raw", 1e7); if (length(b) == 0) break; writeBin(b, out) }
  close(con); close(out)
  md5_of(tmp)
}
paper_files <- function() {
  f <- c(list.files("04_output", recursive = TRUE, full.names = TRUE, all.files = TRUE),
         list.files("05_manuscript", recursive = TRUE, full.names = TRUE, all.files = TRUE))
  f <- f[!grepl("(^|/)not_adopted(/|$)", f) & !grepl("/cache/", f)]
  sort(c(f, "06_diagnostics/fastq_read_names.csv", "06_diagnostics/table1_group_means.csv"), method = "radix")
}
if (length(args) == 2 && args[1] == "base") {
  f <- system2("git", c("ls-files", "04_output", "05_manuscript", "06_diagnostics"), stdout = TRUE)
  m <- data.table(path = f, bytes = file.size(f), md5 = md5_of(f), md5_content = vapply(f, md5_inner, ""))
  fwrite(m, args[2], eol = "\n")
  cat(nrow(m), "files in", args[2], "\n")
} else if (length(args) == 2 && args[1] == "snapshot") {
  f <- paper_files()
  m <- data.table(path = f, bytes = file.size(f), md5 = md5_of(f),
                  mtime = format(file.mtime(f), "%Y-%m-%d %H:%M:%OS6"))
  dirs <- unique(dirname(f))
  m <- rbind(m, data.table(path = paste0(dirs, "/"), bytes = NA_real_, md5 = "folder", mtime = ""))
  fwrite(m, args[2], eol = "\n")
  cat(nrow(m), "files and folders in", args[2], "\n")
} else if (length(args) == 3 && args[1] == "compare") {
  a <- fread(args[2], colClasses = "character"); b <- fread(args[3], colClasses = "character")
  added <- setdiff(b$path, a$path); removed <- setdiff(a$path, b$path)
  both <- merge(a, b, by = "path", suffixes = c("_before", "_after"))
  changed <- both[md5_before != md5_after | mtime_before != mtime_after]
  cat(sprintf("T2: %d paths before, %d after; %d added, %d removed, %d changed (md5 or modification time)\n",
              nrow(a), nrow(b), length(added), length(removed), nrow(changed)))
  if (length(added)) cat("added:\n ", paste(added, collapse = "\n  "), "\n")
  if (length(removed)) cat("removed:\n ", paste(removed, collapse = "\n  "), "\n")
  if (nrow(changed)) print(changed[, .(path, md5_before, md5_after, mtime_before, mtime_after)])
  quit(status = if (length(added) + length(removed) + nrow(changed) == 0) 0 else 1)
} else {
  stop("Usage: Rscript tests/snapshot_paper.R base|snapshot <out.csv>  or  compare <before.csv> <after.csv>", call. = FALSE)
}
