# Which GO.db release reproduces the three "-new" BP tables of the archive?
suppressPackageStartupMessages(library(data.table))
scr <- getwd(); source(file.path(scr, "parse_archive.R"))
A <- "/home/user/WGBS-gigas-ploidy-desiccation/08_archive/GOterms"
newf <- c(PLOIDY = "BOTH-PLOIDY/ploidy-BP-FisherTestResults-new.csv", HEAT = "BOTH-HEAT/HEAT-BP-FisherTestResults-new.csv",
          COMBO = "BOTH-COMBO/COMBO-BP-FisherTestResults-new.csv")
tags <- commandArgs(trailingOnly = TRUE)
res <- list()
for (tag in tags) for (s in names(newf)) {
  fr <- file.path(scr, "runs", tag, sprintf("%s-BP-FisherTestResults.csv", s)); if (!file.exists(fr)) next
  a <- as.data.table(parse_archive(file.path(A, newf[[s]]))); n <- as.data.table(parse_archive(fr))
  m <- merge(a, n, by = "go_id", suffixes = c("_a", "_n"))
  res[[length(res) + 1]] <- data.table(run = tag, set = s, terms_new_table = nrow(a), terms_run = nrow(n),
    only_table = length(setdiff(a$go_id, n$go_id)), only_run = length(setdiff(n$go_id, a$go_id)),
    annotated_differs = sum(m$annotated_a != m$annotated_n), significant_differs = sum(m$significant_a != m$significant_n),
    term_name_differs = sum(m$term_a != m$term_n), printed_p_differs = sum(gsub(" ", "", m$classic_a) != gsub(" ", "", m$classic_n)),
    same_row_order = identical(a$go_id, n$go_id), file_identical = identical(readLines(file.path(A, newf[[s]])), readLines(fr)))
}
options(width = 250); print(rbindlist(res))
