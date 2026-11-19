suppressPackageStartupMessages({library(data.table)})
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
tdir <- file.path(repo, "04_output/05_differential_methylation/tables")
for (md in c("legacy", "corrected")) {
  p0 <- fread(file.path(tdir, sprintf("diffmeth_ploidy_%s.csv.gz", md)), select = c("qvalue", "meth.diff", "united_row"))
  h0 <- fread(file.path(tdir, sprintf("diffmeth_heat_%s.csv.gz", md)), select = c("qvalue", "meth.diff", "united_row"))
  sw <- fread(file.path(tdir, sprintf("dml_counts_by_threshold_%s.csv", md)))
  rec <- rbindlist(lapply(seq(5, 50, 5), function(th) {
    p <- p0[qvalue < 0.01 & abs(meth.diff) > th]; h <- h0[qvalue < 0.01 & abs(meth.diff) > th]
    s <- h[united_row %in% p$united_row]
    data.table(difference = th, set = c("ploidy", "heat", "shared"), n2 = c(nrow(p), nrow(h), nrow(s)),
               n_hyper2 = c(sum(p$meth.diff > 0), sum(h$meth.diff > 0), sum(s$meth.diff > 0)),
               n_hypo2 = c(sum(p$meth.diff < 0), sum(h$meth.diff < 0), sum(s$meth.diff < 0)))
  }))
  m <- merge(sw, rec, by = c("difference", "set"), sort = FALSE)
  cat(md, ": rows", nrow(m), " mismatching rows:", sum(m$n != m$n2 | m$n_hyper != m$n_hyper2 | m$n_hypo != m$n_hypo2), "\n")
  print(m[n != n2 | n_hyper != n_hyper2 | n_hypo != n_hypo2])
  cat("classes sw:", sapply(sw, class), "\n")
}
