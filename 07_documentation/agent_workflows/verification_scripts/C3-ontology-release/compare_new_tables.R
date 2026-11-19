suppressPackageStartupMessages(library(data.table))
scr <- getwd(); source(file.path(scr, "parse_archive.R"))
A <- "/home/user/WGBS-gigas-ploidy-desiccation/08_archive/GOterms"
pairs <- list(PLOIDY = c("BOTH-PLOIDY/PLOIDY-BP-FisherTestResults.csv", "BOTH-PLOIDY/ploidy-BP-FisherTestResults-new.csv"),
              HEAT   = c("BOTH-HEAT/HEAT-BP-FisherTestResults.csv", "BOTH-HEAT/HEAT-BP-FisherTestResults-new.csv"),
              COMBO  = c("BOTH-COMBO/COMBO-BP-FisherTestResults.csv", "BOTH-COMBO/COMBO-BP-FisherTestResults-new.csv"))
for (s in names(pairs)) {
  o <- as.data.table(parse_archive(file.path(A, pairs[[s]][1]))); n <- as.data.table(parse_archive(file.path(A, pairs[[s]][2])))
  m <- merge(o, n, by = "go_id", suffixes = c("_old", "_new"))
  cat(sprintf("\n%s BP: old %d terms, new %d; only old %d, only new %d; shared %d; annotated differs %d; significant differs %d; top10 in common %d; top10 same order %s\n",
              s, nrow(o), nrow(n), length(setdiff(o$go_id, n$go_id)), length(setdiff(n$go_id, o$go_id)), nrow(m),
              sum(m$annotated_old != m$annotated_new), sum(m$significant_old != m$significant_new),
              length(intersect(o$go_id[1:10], n$go_id[1:10])), identical(o$go_id[1:10], n$go_id[1:10])))
  print(cbind(old = paste(o$go_id[1:10], substr(o$term[1:10], 1, 34), o$annotated[1:10], o$significant[1:10], o$classic[1:10]),
              new = paste(n$go_id[1:10], substr(n$term[1:10], 1, 34), n$annotated[1:10], n$significant[1:10], n$classic[1:10])), quote = FALSE)
  cat("p<0.01: old", sum(classic_numeric(o$classic) < 0.01 | grepl("<", o$classic), na.rm = TRUE), " new", sum(classic_numeric(n$classic) < 0.01 | grepl("<", n$classic), na.rm = TRUE), "\n")
}
