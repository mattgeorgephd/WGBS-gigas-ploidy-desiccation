suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv")); setorder(lib, legacy_sample_index)
lc <- fread(file.path(repo, "04_output/04_cpg_coverage/tables/library_coverage.csv"))
fl <- fread(file.path(repo, "04_output/05_differential_methylation/tables/library_filtering_legacy.csv"))
fc <- fread(file.path(repo, "04_output/05_differential_methylation/tables/library_filtering_corrected.csv"))
wq <- function(v, w, p) {  # type-7 quantile of a weighted (count) distribution, exact
  o <- order(v); v <- v[o]; w <- w[o]; n <- sum(w); h <- (n - 1) * p + 1; cw <- cumsum(w)
  lo <- floor(h); hi <- ceiling(h); xlo <- v[which(cw >= lo)[1]]; xhi <- v[which(cw >= hi)[1]]
  xlo + (h - lo) * (xhi - xlo)
}
wmed <- function(v, w) wq(v, w, 0.5)
res <- rbindlist(lapply(seq_len(nrow(lib)), function(j) {
  x <- readLines(sprintf("hist/%s.txt", lib$cov_file[j]))
  h <- fread(text = x[grepl("^hist", x)], header = FALSE, col.names = c("k", "cov", "n"))
  first <- strsplit(sub("^first\t", "", x[grepl("^first", x)]), "\t")[[1]]
  fcov <- as.numeric(first[5]) + as.numeric(first[6])
  nl <- as.numeric(sub("nlines\t", "", x[grepl("^nlines", x)]))
  ge2 <- h[cov >= 2]; cutoff <- wq(ge2$cov, ge2$n, 0.999)
  ge10 <- h[cov >= 10]
  # legacy: first line lost
  hl <- copy(h); hl[cov == fcov, n := n - 1]
  data.table(library_id = lib$library_id[j], nlines = nl, min_cov = min(h$cov), n1 = sum(h$n), n10 = sum(ge10$n), mean_cov = sum(h$cov * h$n) / sum(h$n),
             med_cov = wmed(h$cov, h$n), cutoff = cutoff, n_ge_cutoff = sum(h[cov >= cutoff]$n),
             first_pos = paste0(first[1], ":", first[2]), first_cov = fcov,
             imp_legacy = sum(hl[cov >= 2]$n), imp_corr = sum(ge2$n), n10_legacy = sum(hl[cov >= 10]$n),
             med10 = wmed(ge10$cov, ge10$n), max10 = max(ge10$cov),
             kept_corr = sum(h[cov >= 10 & cov < cutoff]$n), med10_corr = wmed(h[cov >= 10 & cov < cutoff]$cov, h[cov >= 10 & cov < cutoff]$n),
             max_kept_corr = max(h[cov >= 10 & cov < cutoff]$cov),
             cutoff_ge1 = wq(h$cov, h$n, 0.999), cutoff_ge10 = wq(ge10$cov, ge10$n, 0.999))
}))
options(width = 250)
print(res[, .(library_id, nlines, min_cov, n1, n10, mean_cov, med_cov, cutoff, n_ge_cutoff, first_pos, first_cov)])
ok <- function(label, cond) cat(sprintf("%-90s %s\n", label, if (isTRUE(cond)) "OK" else "*** FAIL ***"))
ok("04 library_coverage: order", identical(lc$library_id, res$library_id))
ok("04 n_cpg_1x == file lines", all(lc$n_cpg_1x == res$nlines) && all(res$n1 == res$nlines))
ok("04 min coverage >= 1 in every file", all(res$min_cov >= 1))
ok("04 n_cpg_10x", all(lc$n_cpg_10x == res$n10))
ok("04 mean_coverage", isTRUE(all.equal(lc$mean_coverage, res$mean_cov, tolerance = 1e-12)))
ok("04 median_coverage", all(lc$median_coverage == res$med_cov))
ok("04 upper_cutoff (type-7 quantile among coverage >= 2)", isTRUE(all.equal(lc$upper_cutoff, res$cutoff, tolerance = 1e-12)))
ok("04 n_cpg_at_or_above_cutoff", all(lc$n_cpg_at_or_above_cutoff == res$n_ge_cutoff))
ok("04 first_line_position", identical(lc$first_line_position, res$first_pos))
ok("04 first_line_coverage", all(lc$first_line_coverage == res$first_cov))
ok("05 legacy n_cpg_imported (first line lost)", all(fl$n_cpg_imported == res$imp_legacy))
ok("05 legacy n_cpg_min_coverage", all(fl$n_cpg_min_coverage == res$n10_legacy))
ok("05 legacy n_cpg_kept == n_cpg_min_coverage", all(fl$n_cpg_kept == fl$n_cpg_min_coverage))
ok("05 legacy median_coverage", all(fl$median_coverage == res$med10))
ok("05 legacy max_coverage_kept", all(fl$max_coverage_kept == res$max10))
ok("05 legacy upper_cutoff", isTRUE(all.equal(fl$upper_cutoff, res$cutoff, tolerance = 1e-12)))
ok("05 legacy scaling factor", isTRUE(all.equal(fl$scaling_factor, max(res$med10) / res$med10, tolerance = 1e-12)))
ok("05 corrected n_cpg_imported", all(fc$n_cpg_imported == res$imp_corr))
ok("05 corrected n_cpg_min_coverage", all(fc$n_cpg_min_coverage == res$n10))
ok("05 corrected n_cpg_kept", all(fc$n_cpg_kept == res$kept_corr))
ok("05 corrected n_removed_by_upper", all(fc$n_removed_by_upper == res$n10 - res$kept_corr))
ok("05 corrected median_coverage", all(fc$median_coverage == res$med10_corr))
ok("05 corrected max_coverage_kept", all(fc$max_coverage_kept == res$max_kept_corr))
ok("05 corrected scaling factor", isTRUE(all.equal(fc$scaling_factor, max(res$med10_corr) / res$med10_corr, tolerance = 1e-12)))
cat("\ncutoff if percentile were taken among coverage >= 1 (the legacy coverage script's mincov) and among >= 10:\n")
print(res[, .(library_id, cutoff_ge2 = cutoff, cutoff_ge1, cutoff_ge10)])
fwrite(res, "hist_summary.csv")
