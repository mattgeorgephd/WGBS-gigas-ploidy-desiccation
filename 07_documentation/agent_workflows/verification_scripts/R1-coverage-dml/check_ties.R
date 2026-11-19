suppressPackageStartupMessages({library(data.table)})
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
tdir <- file.path(repo, "04_output/05_differential_methylation/tables")
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv")); setorder(lib, legacy_sample_index)
tp <- as.integer(lib$ploidy == "triploid"); th <- as.integer(lib$treatment == "heat")
for (md in c("legacy", "corrected")) {
  u <- fread(file.path(tdir, sprintf("united_cpg_%s.csv.gz", md)))
  nc <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]); nt <- as.matrix(u[, paste0("numTs", 1:22), with = FALSE])
  for (contrast in c("ploidy", "heat")) {
    d <- fread(file.path(tdir, sprintf("diffmeth_%s_%s.csv.gz", contrast, md)))
    tr <- if (contrast == "ploidy") tp else th
    # recompute exactly as logReg does (tapply sums; meths[2]-meths[1]; 100*)
    ms1 <- rowSums(nc[, tr == 1]); ws1 <- rowSums((nc + nt)[, tr == 1]); ms0 <- rowSums(nc[, tr == 0]); ws0 <- rowSums((nc + nt)[, tr == 0])
    mdiff <- 100 * (ms1 / ws1 - ms0 / ws0)
    for (thr in seq(5, 50, 5)) {
      tie <- which(d$qvalue < 0.01 & abs(d$meth.diff) == thr)
      if (length(tie)) {
        cat(sprintf("[%s %s] threshold %d: %d significant CpG(s) with |meth.diff| written as exactly %d in the CSV\n", md, contrast, thr, length(tie), thr))
        print(data.table(d[tie, .(chr, start, qvalue, meth.diff_csv = sprintf("%.17g", meth.diff))], recomputed = sprintf("%.17g", mdiff[tie]),
                         ms1 = ms1[tie], ws1 = ws1[tie], ms0 = ms0[tie], ws0 = ws0[tie], in_memory_passes = abs(mdiff[tie]) > thr))
      }
    }
    cat(sprintf("[%s %s] rows where sign of (|csv|>thr) differs from (|recomputed|>thr) for any thr in 5..50 (all q): %d\n", md, contrast,
                sum(sapply(seq(5, 50, 5), function(thr) sum((abs(d$meth.diff) > thr) != (abs(mdiff) > thr))))))
    cat(sprintf("[%s %s] rows with meth.diff == 0 in csv: %d ; of those q<0.01: %d\n", md, contrast, sum(d$meth.diff == 0), sum(d$meth.diff == 0 & d$qvalue < 0.01)))
  }
}
