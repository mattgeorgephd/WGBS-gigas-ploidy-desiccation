suppressPackageStartupMessages({library(data.table); library(methylKit)})
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
t5 <- function(f) file.path(repo, "04_output/05_differential_methylation/tables", f)
counts <- fread(t5("dml_counts_by_mode.csv"))
print(t(counts))
ns <- asNamespace("methylKit")
res <- list(); dm <- list(); un <- list()
for (mode in c("legacy", "corrected")) {
  cat("\n================ mode:", mode, "================\n")
  dml <- list(ploidy = fread(t5(sprintf("dml_ploidy_%s.csv", mode))), heat = fread(t5(sprintf("dml_heat_%s.csv", mode))))
  shared <- fread(t5(sprintf("dml_shared_%s.csv", mode)))
  for (contrast in c("ploidy", "heat")) {
    d <- fread(t5(sprintf("diffmeth_%s_%s.csv.gz", contrast, mode)))
    dm[[paste(contrast, mode)]] <- d
    mine <- d[qvalue < 0.01 & abs(meth.diff) > 20]
    cat(sprintf("[%s %s] rows in diffmeth: %d; DML recomputed from q and diff: %d (hyper %d, hypo %d); rows in dml file: %d (hyper %d, hypo %d); same loci: %s\n",
                contrast, mode, nrow(d), nrow(mine), sum(mine$meth.diff > 0), sum(mine$meth.diff < 0), nrow(dml[[contrast]]),
                sum(dml[[contrast]]$direction == "hyper"), sum(dml[[contrast]]$direction == "hypo"),
                identical(paste(mine$chr, mine$start), paste(dml[[contrast]]$chr, dml[[contrast]]$start))))
    cat(sprintf("   direction column consistent with sign: %s; any |diff| exactly 20: %d; NA p: %d; duplicated positions: %d\n",
                all((dml[[contrast]]$meth.diff > 0) == (dml[[contrast]]$direction == "hyper")), sum(abs(d$meth.diff) == 20), sum(is.na(d$pvalue)),
                sum(duplicated(d[, .(chr, start)]))))
    # SLIM q-values recomputed from the stored p-values with methylKit's own functions
    q2 <- ns$QValuesfun(d$pvalue, ns$SLIMfunc(d$pvalue)$pi0_Est)
    cat(sprintf("   SLIM q recomputed from stored p: max abs diff %.3g, max rel diff %.3g; DML with recomputed q: %d\n",
                max(abs(q2 - d$qvalue)), max(abs(q2 - d$qvalue) / pmax(d$qvalue, 1e-300)), sum(q2 < 0.01 & abs(d$meth.diff) > 20)))
    cat(sprintf("   with Benjamini-Hochberg instead of SLIM: %d DML\n", sum(p.adjust(d$pvalue, "BH") < 0.01 & abs(d$meth.diff) > 20)))
  }
  sh <- merge(dml$ploidy[, .(chr, start, dp = meth.diff)], dml$heat[, .(chr, start, dh = meth.diff)], by = c("chr", "start"))
  cls <- sh[, fifelse(dp > 0 & dh > 0, "both_hyper", fifelse(dp < 0 & dh < 0, "both_hypo", "opposite"))]
  cat(sprintf("[shared %s] recomputed overlap: %d (both hyper %d, both hypo %d, opposite %d); file: %d (%s)\n", mode, nrow(sh),
              sum(cls == "both_hyper"), sum(cls == "both_hypo"), sum(cls == "opposite"), nrow(shared),
              paste(names(table(shared$class)), table(shared$class), collapse = ", ")))
  cat("   same loci as file:", setequal(paste(sh$chr, sh$start), paste(shared$chr, shared$start)), "\n")
  mode_i <- mode; row <- counts[counts$mode == mode_i]
  stopifnot(nrow(row) == 1)
  cat(sprintf("   counts file row agrees: cpg %s, ploidy %s, heat %s, shared %s\n",
              row$cpg_analysed == nrow(dm[[paste("ploidy", mode)]]),
              all(c(row$ploidy_dml, row$ploidy_hyper, row$ploidy_hypo) == c(nrow(dml$ploidy), sum(dml$ploidy$meth.diff > 0), sum(dml$ploidy$meth.diff < 0))),
              all(c(row$heat_dml, row$heat_hyper, row$heat_hypo) == c(nrow(dml$heat), sum(dml$heat$meth.diff > 0), sum(dml$heat$meth.diff < 0))),
              all(c(row$shared, row$shared_both_hyper, row$shared_both_hypo, row$shared_opposite) ==
                    c(nrow(sh), sum(cls == "both_hyper"), sum(cls == "both_hypo"), sum(cls == "opposite")))))
  # regions: chain DML of one list whose neighbours are <= 1000 bp apart
  runs <- function(d) { d <- d[order(chr, start)]; g <- cumsum(c(TRUE, d$chr[-1] != d$chr[-nrow(d)] | diff(d$start) > 1000)); tab <- table(g)
    c(n_regions = length(tab), largest = max(tab), singles = sum(tab == 1), in5 = sum(tab[tab >= 5]),
      pct_with_neighbour_within_1kb = round(100 * mean(tab[as.character(g)] > 1), 1)) }
  cat("   regions ploidy:", paste(names(runs(dml$ploidy)), runs(dml$ploidy), collapse = "; "), "\n")
  cat("   regions heat:  ", paste(names(runs(dml$heat)), runs(dml$heat), collapse = "; "), "\n")
  cat("   regions shared:", paste(names(runs(shared)), runs(shared), collapse = "; "), "\n")
  # the locus
  inloc <- function(d) d[chr == "NC_047563.1" & start >= 39633323 - 1 & start <= 39635520, .N]
  cat(sprintf("   locus NC_047563.1:39,633,323-39,635,520: tested CpGs %d, ploidy DML %d, heat DML %d, shared %d\n",
              inloc(dm[[paste("ploidy", mode)]]), inloc(dml$ploidy), inloc(dml$heat), inloc(shared)))
  res[[mode]] <- list(dml = dml, shared = shared)
}
# ---- legacy against corrected ----
cat("\n================ legacy against corrected ================\n")
for (contrast in c("ploidy", "heat")) {
  L <- dm[[paste(contrast, "legacy")]]; C <- dm[[paste(contrast, "corrected")]]
  m <- merge(L, C, by = c("chr", "start"), suffixes = c("_L", "_C"))
  cat(sprintf("[%s] legacy CpGs %d, corrected %d, in both %d; corrected not in legacy %d; legacy not in corrected %d\n", contrast,
              nrow(L), nrow(C), nrow(m), nrow(C) - nrow(m), nrow(L) - nrow(m)))
  cat(sprintf("   common CpGs: p-values identical %s (max abs diff %.3g); meth.diff identical %s (max abs diff %.3g); q ratio corrected/legacy range %.4f to %.4f\n",
              identical(m$pvalue_L, m$pvalue_C), max(abs(m$pvalue_L - m$pvalue_C)), identical(m$meth.diff_L, m$meth.diff_C),
              max(abs(m$meth.diff_L - m$meth.diff_C)), min(m$qvalue_C / m$qvalue_L), max(m$qvalue_C / m$qvalue_L)))
  dl <- res$legacy$dml[[contrast]]; dc <- res$corrected$dml[[contrast]]
  kl <- paste(dl$chr, dl$start); kc <- paste(dc$chr, dc$start); kC <- paste(C$chr, C$start)
  cat(sprintf("   legacy DML %d: %d still DML, %d removed with their CpG by the filter, %d still tested but no longer DML; corrected DML not in legacy list: %d\n",
              length(kl), sum(kl %in% kc), sum(!kl %in% kC), sum(kl %in% kC & !kl %in% kc), sum(!kc %in% kl)))
}
sl <- res$legacy$shared; sc <- res$corrected$shared
cat(sprintf("[shared] legacy %d, corrected %d; lost %d; gained %d; opposite lost: %d of %d (%.1f%%)\n", nrow(sl), nrow(sc),
            sum(!paste(sl$chr, sl$start) %in% paste(sc$chr, sc$start)), sum(!paste(sc$chr, sc$start) %in% paste(sl$chr, sl$start)),
            sum(sl$class == "opposite") - sum(sc$class == "opposite"), sum(sl$class == "opposite"),
            100 * (1 - sum(sc$class == "opposite") / sum(sl$class == "opposite"))))
# where the lost DML are: regions of removed legacy DML
for (contrast in c("ploidy", "heat")) {
  dl <- res$legacy$dml[[contrast]]; C <- dm[[paste(contrast, "corrected")]]
  lost <- dl[!paste(chr, start) %in% paste(C$chr, C$start)][order(chr, start)]
  lost[, g := cumsum(c(TRUE, chr[-1] != chr[-.N] | diff(start) > 1000))]
  top <- lost[, .(chr = chr[1], from = min(start), to = max(start), n = .N), by = g][order(-n)][1:6]
  cat(sprintf("[%s] legacy DML at CpGs removed by the filter: %d, in %d runs; largest runs:\n", contrast, nrow(lost), uniqueN(lost$g)))
  print(top[, .(chr, from, to, n)])
  cat("   by sequence class:", paste(names(table(substr(lost$chr, 1, 2))), table(substr(lost$chr, 1, 2)), collapse = ", "),
      "; mitochondrial:", sum(lost$chr == "NC_001276.1"), "\n")
}
# CpGs removed by the filter, by sequence
L <- dm[["ploidy legacy"]]; C <- dm[["ploidy corrected"]]
rem <- L[!paste(chr, start) %in% paste(C$chr, C$start)]
cat(sprintf("\nCpGs removed by the filter: %d (%.2f%% of %d); mitochondrial %d; unplaced scaffolds %d; linkage groups %d\n", nrow(rem),
            100 * nrow(rem) / nrow(L), nrow(L), sum(rem$chr == "NC_001276.1"), sum(grepl("^NW_", rem$chr)), sum(grepl("^NC_0475", rem$chr))))
cat("tested CpGs on the mitochondrial genome: legacy", sum(L$chr == "NC_001276.1"), " corrected", sum(C$chr == "NC_001276.1"), "\n")
saveRDS(list(res = res), "dml_check.rds")
