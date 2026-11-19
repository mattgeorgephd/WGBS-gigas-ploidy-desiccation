# Date the earlier archived topGO tables by their term sets and Annotated counts,
# which depend only on the universe file and the GO.db release.
suppressPackageStartupMessages(library(data.table))
scr <- getwd(); source(file.path(scr, "parse_archive.R"))
R <- "/home/user/WGBS-gigas-ploidy-desiccation/08_archive"
# archived table, universe key, ontology, gannet date of the table
tabs <- rbind(
  data.table(file = "tables/ploidy-CC-FisherTestResults.csv", universe = "july", ont = "CC", dated = "2021-07-29"),
  data.table(file = "tables/ploidy-MF-FisherTestResults.csv", universe = "july", ont = "MF", dated = "2021-07-29"),
  data.table(file = "GOterms/BOTH-COMBO/ploidy-BP-FisherTestResults.csv", universe = c("july", "ronit"), ont = "BP", dated = "2021-11-21"),
  data.table(file = "GOterms/BOTH-COMBO/ploidy-CC-FisherTestResults.csv", universe = c("july", "ronit"), ont = "CC", dated = "2021-11-21"),
  data.table(file = "GOterms/BOTH-COMBO/ploidy-MF-FisherTestResults.csv", universe = c("july", "ronit"), ont = "MF", dated = "2021-11-21"),
  data.table(file = "GOterms/RONIT-PLOIDY/PLOIDY-BP-FisherTestResults.csv", universe = "ronit", ont = "BP", dated = "2021-11-29"),
  data.table(file = "GOterms/RONIT-PLOIDY/PLOIDY-CC-FisherTestResults.csv", universe = "ronit", ont = "CC", dated = "2021-11-29"),
  data.table(file = "GOterms/RONIT-PLOIDY/PLOIDY-MF-FisherTestResults.csv", universe = "ronit", ont = "MF", dated = "2021-11-29"),
  data.table(file = "GOterms/RONIT-HEAT/HEAT-BP-FisherTestResults.csv", universe = "ronit", ont = "BP", dated = "2021-11-29"),
  data.table(file = "GOterms/RONIT-COMBO/COMBO-MF-FisherTestResults.csv", universe = "ronit", ont = "MF", dated = "2021-11-29"),
  data.table(file = "GOterms/YAAMINI-PLOIDY/ploidy-CC-FisherTestResults.csv", universe = "yaamini", ont = "CC", dated = "2021-11-29"),
  data.table(file = "GOterms/YAAMINI-PLOIDY/ploidy-MF-FisherTestResults.csv", universe = "yaamini", ont = "MF", dated = "2021-11-29"))
rel <- c(godb3114 = "3.11.4", godb3121 = "3.12.1", godb313 = "3.13.0", godb314 = "3.14.0")
out <- list()
for (i in seq_len(nrow(tabs))) {
  l <- readLines(file.path(R, tabs$file[i])); l <- sub(",$", "", l)          # one hand-edited file has a trailing comma
  tf <- tempfile(tmpdir = file.path(scr, "tmp")); writeLines(l, tf); a <- as.data.table(parse_archive(tf)); unlink(tf)
  for (r in names(rel)) {
    f <- file.path(scr, "universe_terms", r, sprintf("%s_%s.csv", tabs$universe[i], tabs$ont[i])); if (!file.exists(f)) next
    u <- fread(f)
    m <- merge(a[, .(go_id, annotated)], u, by = "go_id", suffixes = c("_table", "_universe"))
    out[[length(out) + 1]] <- data.table(table = tabs$file[i], dated = tabs$dated[i], universe = tabs$universe[i], GO.db = rel[[r]],
      terms_table = nrow(a), terms_universe = nrow(u), table_ids_absent = length(setdiff(a$go_id, u$go_id)),
      extra_ids = length(setdiff(u$go_id, a$go_id)), annotated_differs = sum(m$annotated_table != m$annotated_universe),
      exact = setequal(a$go_id, u$go_id) && all(m$annotated_table == m$annotated_universe))
  }
}
res <- rbindlist(out); options(width = 250)
print(res)
fwrite(res, file.path(scr, "universe_dating_results.csv"))
cat("\nReleases that reproduce each table exactly (ids and Annotated):\n")
print(res[, .(exact_under = paste(GO.db[exact], collapse = ", ")), by = .(table, dated, universe)])
