suppressPackageStartupMessages({library(data.table)})
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
tdir <- file.path(repo, "04_output/05_differential_methylation/tables")
bdir <- file.path(repo, "04_output/05_differential_methylation/bed")
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv"))
setorder(lib, legacy_sample_index)
tp <- as.integer(lib$ploidy == "triploid"); th <- as.integer(lib$treatment == "heat")
ok <- function(label, cond) cat(sprintf("%-90s %s\n", label, if (isTRUE(cond)) "OK" else "*** FAIL ***"))
for (mode in c("legacy", "corrected")) {
  cat("\n=========== mode:", mode, "===========\n")
  u <- fread(file.path(tdir, sprintf("united_cpg_%s.csv.gz", mode)))
  cat("united rows:", nrow(u), " cols:", ncol(u), "\n")
  print(sapply(u[, 1:10], class))
  ok("united: 4 + 66 columns", ncol(u) == 70)
  ok("united: end == start + 2 everywhere", all(u$end == u$start + 2L))
  ok("united: strand all '*'", all(u$strand == "*"))
  ok("united: no duplicated (chr,start)", !anyDuplicated(u[, .(chr, start)]))
  # sorted by chr (C locale) then start?
  o <- order(u$chr, u$start, method = "radix")
  ok("united: sorted by chr (C-locale radix) then start", identical(o, seq_len(nrow(u))))
  covc <- paste0("coverage", 1:22); cc <- paste0("numCs", 1:22); tc <- paste0("numTs", 1:22)
  ok("united: column names coverage/numCs/numTs 1..22 in order",
     identical(names(u)[-(1:4)], as.vector(rbind(covc, cc, tc))))
  cv <- as.matrix(u[, ..covc]); nc <- as.matrix(u[, ..cc]); nt <- as.matrix(u[, ..tc])
  cat("  rows where numCs+numTs != coverage (rounding), by library:\n   "); cat(colSums(nc + nt != cv), "\n")
  cat("  max |numCs+numTs-coverage|:", max(abs(nc + nt - cv)), "\n")
  cat("  min coverage per library:\n   "); cat(apply(cv, 2, min), "\n")
  cat("  max coverage per library:\n   "); cat(apply(cv, 2, max), "\n")
  filt <- fread(file.path(tdir, sprintf("library_filtering_%s.csv", mode)))
  ok("filtering rows in sample-sheet order", identical(filt$library_id, lib$library_id))
  # minimum scaled coverage must be round(sf*10)
  ok("min scaled coverage == round(sf*10) for every library", all(apply(cv, 2, min) == round(filt$scaling_factor * 10)))
  if (mode == "corrected") {
    lim <- round(filt$scaling_factor * (ceiling(filt$upper_cutoff) - 1))
    cat("  corrected: max allowed scaled coverage per library:\n   "); cat(lim, "\n")
    ok("corrected: max scaled coverage <= round(sf*(ceil(cutoff)-1))", all(apply(cv, 2, max) <= lim))
  }
  # wmean differences
  wdiff <- function(tr) 100 * (rowSums(nc[, tr == 1]) / rowSums((nc + nt)[, tr == 1]) - rowSums(nc[, tr == 0]) / rowSums((nc + nt)[, tr == 0]))
  dm <- list()
  for (contrast in c("ploidy", "heat")) {
    d <- fread(file.path(tdir, sprintf("diffmeth_%s_%s.csv.gz", contrast, mode)))
    ok(sprintf("diffmeth_%s: same rows/order as united", contrast), identical(d$chr, u$chr) && identical(d$start, u$start) && identical(d$end, u$end))
    ok(sprintf("diffmeth_%s: united_row == 1..n", contrast), identical(d$united_row, seq_len(nrow(u))))
    tr <- if (contrast == "ploidy") tp else th
    md <- wdiff(tr)
    cat(sprintf("  diffmeth_%s: max |meth.diff - recomputed (group1 - group0, pooled)| = %.3g\n", contrast, max(abs(d$meth.diff - md))))
    ok(sprintf("diffmeth_%s: pvalue in [0,1], qvalue finite", contrast), all(d$pvalue >= 0 & d$pvalue <= 1) && all(is.finite(d$qvalue)))
    cat(sprintf("  diffmeth_%s: NA p: %d ; p==1: %d ; q range %.3g..%.3g\n", contrast, sum(is.na(d$pvalue)), sum(d$pvalue == 1), min(d$qvalue), max(d$qvalue)))
    # monotone q in p
    op <- order(d$pvalue)
    ok(sprintf("diffmeth_%s: q non-decreasing in p", contrast), !is.unsorted(d$qvalue[op]))
    dml <- fread(file.path(tdir, sprintf("dml_%s_%s.csv", contrast, mode)))
    exp <- d[qvalue < 0.01 & abs(meth.diff) > 20]
    ok(sprintf("dml_%s: same loci as diffmeth filter (n=%d vs %d)", contrast, nrow(dml), nrow(exp)),
       nrow(dml) == nrow(exp) && identical(dml$united_row, exp$united_row) && identical(dml$chr, exp$chr) && identical(dml$start, exp$start))
    ok(sprintf("dml_%s: meth.diff/q identical to diffmeth rows", contrast), isTRUE(all.equal(dml$meth.diff, exp$meth.diff, tolerance = 1e-12)) && isTRUE(all.equal(dml$qvalue, exp$qvalue, tolerance = 1e-12)))
    ok(sprintf("dml_%s: direction == sign(meth.diff)", contrast), all(dml$direction == ifelse(dml$meth.diff > 0, "hyper", "hypo")))
    cat(sprintf("  dml_%s: n=%d hyper=%d hypo=%d ; min|diff|=%.4f ; max q=%.3g\n", contrast, nrow(dml), sum(dml$direction == "hyper"), sum(dml$direction == "hypo"), min(abs(dml$meth.diff)), max(dml$qvalue)))
    cat(sprintf("  diffmeth_%s: rows with q<0.01 and |diff| in (19.999,20.001): %d ; q in (0.0099,0.0101) & |diff|>20: %d\n", contrast,
        sum(d$qvalue < 0.01 & abs(abs(d$meth.diff) - 20) < 1e-3), sum(abs(d$qvalue - 0.01) < 1e-4 & abs(d$meth.diff) > 20)))
    bed <- fread(file.path(bdir, sprintf("dml_%s_%s.bed", contrast, mode)), header = FALSE)
    ok(sprintf("bed %s: same loci and difference", contrast), nrow(bed) == nrow(dml) && identical(bed$V1, dml$chr) && identical(bed$V2, dml$start) && identical(bed$V3, dml$end) && isTRUE(all.equal(bed$V4, dml$meth.diff, tolerance = 1e-9)))
    dm[[contrast]] <- list(all = d, dml = dml)
  }
  sh <- fread(file.path(tdir, sprintf("dml_shared_%s.csv", mode)))
  ex <- intersect(dm$ploidy$dml$united_row, dm$heat$dml$united_row)
  ok(sprintf("shared: equals intersection of united_row (n=%d vs %d)", nrow(sh), length(ex)), identical(sort(sh$united_row), sort(ex)))
  ok("shared: sorted by united_row", !is.unsorted(sh$united_row))
  ok("shared: meth.diff_ploidy / meth.diff_heat equal diffmeth rows",
     isTRUE(all.equal(sh$meth.diff_ploidy, dm$ploidy$all$meth.diff[sh$united_row], tolerance = 1e-12)) && isTRUE(all.equal(sh$meth.diff_heat, dm$heat$all$meth.diff[sh$united_row], tolerance = 1e-12)))
  ok("shared: chr/start equal united rows", identical(sh$chr, u$chr[sh$united_row]) && identical(sh$start, u$start[sh$united_row]))
  cls <- ifelse(sh$meth.diff_ploidy > 0 & sh$meth.diff_heat > 0, "both_hyper", ifelse(sh$meth.diff_ploidy < 0 & sh$meth.diff_heat < 0, "both_hypo", "opposite"))
  ok("shared: class consistent with signs", identical(cls, sh$class))
  print(table(sh$class))
  cat("  shared: heat hyper", sum(sh$meth.diff_heat > 0), " heat hypo", sum(sh$meth.diff_heat < 0), "\n")
  print(names(sh))
  bed <- fread(file.path(bdir, sprintf("dml_shared_%s.bed", mode)), header = FALSE)
  ok("bed shared: same loci, heat difference", nrow(bed) == nrow(sh) && identical(bed$V1, sh$chr) && identical(bed$V2, sh$start) && identical(bed$V3, sh$end) && isTRUE(all.equal(bed$V4, sh$meth.diff_heat, tolerance = 1e-9)))
  # sweep
  sw <- fread(file.path(tdir, sprintf("dml_counts_by_threshold_%s.csv", mode)))
  rec <- rbindlist(lapply(seq(5, 50, 5), function(th) {
    p <- dm$ploidy$all[qvalue < 0.01 & abs(meth.diff) > th]; h <- dm$heat$all[qvalue < 0.01 & abs(meth.diff) > th]
    s <- h[united_row %in% p$united_row]
    data.table(mode = mode, difference = th, set = c("ploidy", "heat", "shared"), n = c(nrow(p), nrow(h), nrow(s)),
               n_hyper = c(sum(p$meth.diff > 0), sum(h$meth.diff > 0), sum(s$meth.diff > 0)),
               n_hypo = c(sum(p$meth.diff < 0), sum(h$meth.diff < 0), sum(s$meth.diff < 0)))
  }))
  ok("sweep: counts equal recomputation from diffmeth files", identical(sw$set, rec$set) && all(sw$difference == rec$difference) && all(sw$n == rec$n) && all(sw$n_hyper == rec$n_hyper) && all(sw$n_hypo == rec$n_hypo))
  ok("sweep: pct_hyper == 100*n_hyper/n", isTRUE(all.equal(sw$pct_hyper[sw$n > 0], (100 * sw$n_hyper / sw$n)[sw$n > 0])) && all(is.na(sw$pct_hyper[sw$n == 0])))
  ok("sweep: n non-increasing with threshold in each set", all(sw[, all(diff(n) <= 0), by = set]$V1))
  ok("sweep: n_hyper + n_hypo == n", all(sw$n_hyper + sw$n_hypo == sw$n))
  hd <- as.data.frame(fread(file.path(tdir, "dml_counts_by_mode.csv"))); hd <- hd[hd$mode == mode, ]
  ok("headline: equals files", hd$cpg_analysed == nrow(u) && hd$ploidy_dml == nrow(dm$ploidy$dml) && hd$heat_dml == nrow(dm$heat$dml) && hd$shared == nrow(sh) &&
       hd$shared_both_hyper == sum(sh$class == "both_hyper") && hd$shared_both_hypo == sum(sh$class == "both_hypo") && hd$shared_opposite == sum(sh$class == "opposite") &&
       hd$ploidy_hyper == sum(dm$ploidy$dml$meth.diff > 0) && hd$heat_hyper == sum(dm$heat$dml$meth.diff > 0))
  saveRDS(list(pp = dm$ploidy$all$pvalue, ph = dm$heat$all$pvalue, qp = dm$ploidy$all$qvalue, qh = dm$heat$all$qvalue), sprintf("pq_%s.rds", mode))
  rm(u, cv, nc, nt, dm); gc(verbose = FALSE)
}
