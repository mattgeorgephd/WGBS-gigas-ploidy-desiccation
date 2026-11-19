# Compare one re-run (runs/<tag>/) with the archived tables.
# Usage: Rscript compare_run.R <tag>
suppressPackageStartupMessages(library(data.table))
tag <- commandArgs(trailingOnly = TRUE)[1]
scr  <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C3-ontology-release"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
source(file.path(scr, "parse_archive.R"))
run <- file.path(scr, "runs", tag)
cat("== run", tag, "==\n"); cat(readLines(file.path(run, "versions.tsv")), sep = "\n")
res <- list()
for (f in sort(list.files(run, pattern = "FisherTestResults\\.csv$"))) {
  s <- sub("-.*", "", f); o <- sub("^[A-Z]+-([A-Z]+)-.*", "\\1", f)
  fa <- file.path(repo, "08_archive/GOterms", paste0("BOTH-", s), f)
  a <- as.data.table(parse_archive(fa)); n <- as.data.table(parse_archive(file.path(run, f)))
  la <- sub("\r$", "", readLines(fa)); ln <- readLines(file.path(run, f))
  m <- merge(a, n, by = "go_id", suffixes = c("_a", "_n"))
  res[[length(res) + 1]] <- data.table(
    set = s, ont = o, terms_archive = nrow(a), terms_run = nrow(n),
    only_archive = length(setdiff(a$go_id, n$go_id)), only_run = length(setdiff(n$go_id, a$go_id)), shared = nrow(m),
    shared_annotated_differs = sum(m$annotated_a != m$annotated_n),
    shared_significant_differs = sum(m$significant_a != m$significant_n),
    shared_term_name_differs = sum(m$term_a != m$term_n),
    shared_printed_p_differs = sum(gsub(" ", "", m$classic_a) != gsub(" ", "", m$classic_n)),
    top10_ids_in_common = length(intersect(a$go_id[1:10], n$go_id[1:10])),
    top10_identical_order = identical(a$go_id[1:10], n$go_id[1:10]),
    same_row_order = identical(a$go_id, n$go_id),
    file_identical = identical(la, ln),
    lines_differing = if (length(la) == length(ln)) sum(la != ln) else NA_integer_)
}
res <- rbindlist(res)
options(width = 260)
print(res)
fwrite(cbind(tag = tag, res), file.path(run, "comparison_with_archive.csv"))
