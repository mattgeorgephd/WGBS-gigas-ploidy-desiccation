suppressPackageStartupMessages(library(data.table))
scr <- getwd(); source(file.path(scr, "parse_archive.R"))
R <- "/home/user/WGBS-gigas-ploidy-desiccation/08_archive"; A <- file.path(R, "GOterms")
get_revigo <- function(f) {
  l <- readLines(f); l <- l[grepl('c\\("GO:[0-9]{7}"', l)]
  id <- sub('.*c\\("(GO:[0-9]{7})".*', "\\1", l)
  # the "value" column, when present, is the 7th field
  nm <- eval(parse(text = grep("^revigo.names", readLines(f), value = TRUE)[1])); 
  vals <- lapply(l, function(x) eval(parse(text = sub("^revigo.data <- rbind\\(", "", sub("\\)?[,;]?\\)?;?$", ")", sub(",$", "", x))))))
  d <- as.data.table(do.call(rbind, vals)); setnames(d, nm[seq_len(ncol(d))]); d
}
tabs <- list(
  ploidy_old = "BOTH-PLOIDY/PLOIDY-BP-FisherTestResults.csv", ploidy_new = "BOTH-PLOIDY/ploidy-BP-FisherTestResults-new.csv",
  heat_old = "BOTH-HEAT/HEAT-BP-FisherTestResults.csv", heat_new = "BOTH-HEAT/HEAT-BP-FisherTestResults-new.csv",
  combo_old = "BOTH-COMBO/COMBO-BP-FisherTestResults.csv", combo_new = "BOTH-COMBO/COMBO-BP-FisherTestResults-new.csv",
  combo_ploidy_2021 = "BOTH-COMBO/ploidy-BP-FisherTestResults.csv",
  yaamini_ploidy = "YAAMINI-PLOIDY/ploidy-BP-FisherTestResults.csv", ronit_ploidy = "RONIT-PLOIDY/PLOIDY-BP-FisherTestResults.csv",
  ronit_heat = "RONIT-HEAT/HEAT-BP-FisherTestResults.csv", ronit_combo = "RONIT-COMBO/COMBO-BP-FisherTestResults.csv")
tabs <- lapply(tabs, function(f) { d <- as.data.table(parse_archive(file.path(A, f))); d[, p := classic_numeric(classic)]; d[grepl("<", classic), p := 1e-30]; d })
for (f in c("figures/revigo_figure/Revigo-PLOIDY-BP.R", "code/Revigo.R")) {
  rv <- get_revigo(file.path(R, f))
  cat(sprintf("\n%s: %d terms; columns: %s\n", f, nrow(rv), paste(names(rv), collapse = ",")))
  for (n in names(tabs)) {
    t <- tabs[[n]]; sig <- t[p < 0.01]
    line <- sprintf("  %-18s table terms %5d, p<0.01 %4d | revigo ids in table %3d, in p<0.01 set %3d", n, nrow(t), nrow(sig),
                    sum(rv$term_ID %in% t$go_id), sum(rv$term_ID %in% sig$go_id))
    if ("value" %in% names(rv)) {
      m <- merge(rv[, .(go_id = term_ID, value = as.numeric(value))], t[, .(go_id, p)], by = "go_id")
      line <- paste0(line, sprintf(" | log10(p) equal to revigo value: %d of %d", sum(abs(m$value - log10(m$p)) < 1e-6), nrow(m)))
    }
    cat(line, "\n")
  }
}
