# Independent odds ratios and share percentages from my own overlap table.
suppressPackageStartupMessages(library(data.table))
scr  <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C5-feature-odds"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
cpg <- readRDS(file.path(scr, "cpg_features.rds"))
N <- nrow(cpg)

## ---- eligibility -----------------------------------------------------------
cat("== eligibility (mean of per-library percentages) ==\n")
cat("n mean < 10:", sum(cpg$mean_unw < 10), sprintf("(%.2f%%)", 100 * mean(cpg$mean_unw < 10)), "\n")
cat("n mean > 90:", sum(cpg$mean_unw > 90), sprintf("(%.2f%%)", 100 * mean(cpg$mean_unw > 90)), "\n")
for (def in c("unw", "pooled")) {
  m <- cpg[[paste0("mean_", def)]]
  cat(def, ": [10,90] closed:", sum(m >= 10 & m <= 90),
      " (10,90) open:", sum(m > 10 & m < 90),
      " [10,90):", sum(m >= 10 & m < 90), " (10,90]:", sum(m > 10 & m <= 90), "\n")
}
cpg[, elig := mean_unw >= 10 & mean_unw <= 90]
cat("eligible:", sum(cpg$elig), sprintf("(%.2f%%)", 100 * mean(cpg$elig)), "\n")
cat("ploidy DML not eligible:", sum(cpg$dml_ploidy & !cpg$elig),
    " heat DML not eligible:", sum(cpg$dml_heat & !cpg$elig),
    " either-list DML not eligible (unique loci):", sum((cpg$dml_ploidy | cpg$dml_heat) & !cpg$elig), "\n")
print(cpg[(dml_ploidy | dml_heat) & !elig,
          .(chr, start, mean_unw, mean_pooled, dml_ploidy, dml_heat, md_ploidy, md_heat, in_gene, in_CDS, in_intron, in_TE, in_intergenic)])
cat("range of mean_unw among ploidy DML:", range(cpg$mean_unw[cpg$dml_ploidy]), "\n")
cat("range of mean_unw among heat DML:", range(cpg$mean_unw[cpg$dml_heat]), "\n")

## ---- overlap counts per set ------------------------------------------------
tr <- sub("^in_", "", grep("^in_", names(cpg), value = TRUE))
sets <- list(
  "all analysed CpGs" = rep(TRUE, N),
  "eligible CpGs" = cpg$elig,
  "ploidy DML" = cpg$dml_ploidy,
  "heat DML" = cpg$dml_heat,
  "shared" = cpg$dml_shared,
  "shared, both hyper" = cpg$shared_class %in% "both hyper",
  "shared, both hypo" = cpg$shared_class %in% "both hypo",
  "shared, opposite" = cpg$shared_class %in% "opposite"
)
eight <- c("CDS", "exon_UTR", "intron", "lncRNA", "upstream_flank", "downstream_flank", "intergenic", "TE")
ov <- rbindlist(lapply(names(sets), function(s) {
  i <- sets[[s]]
  n <- vapply(tr, function(t) sum(cpg[[paste0("in_", t)]][i]), numeric(1))
  d <- data.table(set = s, n_loci = sum(i), track = tr, n_overlap = n)
  d[, pct_of_loci := 100 * n_overlap / n_loci]
  d[track %in% eight, share8 := 100 * n_overlap / sum(n_overlap)]
  d
}))
fwrite(ov, file.path(scr, "my_overlap_by_track.csv"))

pipe <- fread(file.path(repo, "04_output/07_genome_features/tables/feature_overlap_by_track_legacy.csv"))
cmp <- merge(ov, pipe, by = c("set", "track"), all = TRUE, suffixes = c("_mine", "_pipe"))
cat("\n== comparison with pipeline feature_overlap_by_track_legacy.csv ==\n")
cat("rows mine:", nrow(ov), " rows pipeline:", nrow(pipe), " merged:", nrow(cmp), "\n")
cat("n_loci identical:", all(cmp$n_loci_mine == cmp$n_loci_pipe), "\n")
cat("n_overlap identical:", all(cmp$n_overlap_mine == cmp$n_overlap_pipe), "\n")
cat("max |share8 - legacy_share_pct|:", max(abs(cmp$share8 - cmp$legacy_share_pct), na.rm = TRUE), "\n")
cat("NA pattern identical:", identical(is.na(cmp$share8), is.na(cmp$legacy_share_pct)), "\n")
print(cmp[n_overlap_mine != n_overlap_pipe | n_loci_mine != n_loci_pipe])

## ---- Fisher tests ----------------------------------------------------------
fish <- function(a, b, c, d) {
  # a: DML in feature, b: DML not in feature, c: comparison in feature, d: comparison not in
  ft <- fisher.test(matrix(c(a, c, b, d), 2))
  data.table(a = a, b = b, c = c, d = d,
             or_cmle = unname(ft$estimate), lo = ft$conf.int[1], hi = ft$conf.int[2],
             or_sample = (a * d) / (b * c), p = ft$p.value)
}
variants <- function(dml, bg, feat) {
  # V1: DML (all of the list) vs background CpGs that are not in this DML list
  # V2: DML restricted to the background vs background CpGs not in the list
  # V3: DML (all) vs the whole background including the DML themselves
  other <- bg & !dml
  v1 <- fish(sum(dml & feat), sum(dml & !feat), sum(other & feat), sum(other & !feat))
  dmlb <- dml & bg
  v2 <- fish(sum(dmlb & feat), sum(dmlb & !feat), sum(other & feat), sum(other & !feat))
  v3 <- fish(sum(dml & feat), sum(dml & !feat), sum(bg & feat), sum(bg & !feat))
  rbind(cbind(variant = "V1 all DML vs other background", v1),
        cbind(variant = "V2 DML in background vs other background", v2),
        cbind(variant = "V3 all DML vs whole background", v3))
}
test_tracks <- c("gene", "intron", "CDS", "TE", "intergenic", "exon_UTR", "lncRNA", "upstream_flank", "downstream_flank",
                 "exon", "mRNA", "flanks", "nonCDS")
bgs <- list("all analysed CpGs" = rep(TRUE, N), "eligible CpGs" = cpg$elig)
res <- rbindlist(lapply(names(bgs), function(b) rbindlist(lapply(c("ploidy DML", "heat DML", "shared"), function(s)
  rbindlist(lapply(test_tracks, function(t)
    cbind(background = b, set = s, track = t, variants(sets[[s]], bgs[[b]], cpg[[paste0("in_", t)]]))))))))
fwrite(res, file.path(scr, "my_fisher_variants.csv"))

pa <- fread(file.path(repo, "04_output/07_genome_features/tables/feature_association_tests_legacy.csv"))
cat("\n== which 2x2 variant does the pipeline table match? ==\n")
for (v in unique(res$variant)) {
  m <- merge(res[variant == v], pa, by = c("background", "set", "track"))
  cat(sprintf("%-45s rows %d  max|OR diff| %.3g  max rel diff %.3g  max|pct_other diff| %.3g  dml_total equal %s\n",
              v, nrow(m), max(abs(m$or_cmle - m$odds_ratio)), max(abs(m$or_cmle / m$odds_ratio - 1)),
              max(abs(100 * m$c / (m$c + m$d) - m$pct_other_cpg)), all(m$a + m$b == m$dml_total)))
}

cat("\n== claim table (V1 = DML vs other CpGs of the background) ==\n")
claim <- data.table(
  track = rep(c("gene", "intron", "CDS", "TE", "intergenic"), each = 4),
  set = rep(c("ploidy DML", "ploidy DML", "heat DML", "heat DML"), 5),
  background = rep(c("all analysed CpGs", "eligible CpGs"), 10),
  claimed = c(7.1, 2.5, 6.1, 2.2, 1.9, 0.80, 2.0, 0.87, 1.5, 1.7, 1.4, 1.6, 0.77, 0.54, 0.83, 0.57, 0.11, 0.30, 0.14, 0.37))
for (v in unique(res$variant)) {
  m <- merge(claim, res[variant == v], by = c("background", "set", "track"), sort = FALSE)
  m[, rounded := signif(or_cmle, 2)]
  m[, rounded_sample := signif(or_sample, 2)]
  m[, ok := abs(rounded - claimed) < 1e-9]
  cat("\n", v, " -- agree at 2 significant figures:", sum(m$ok), "of", nrow(m), "\n")
  print(m[, .(track, set, background, claimed, or_cmle = round(or_cmle, 4), or_sample = round(or_sample, 4),
              lo = round(lo, 3), hi = round(hi, 3), a, b, c, d, p = signif(p, 3), ok)])
}

cat("\n== other tracks against eligible (V1) ==\n")
print(res[variant == "V1 all DML vs other background" & track %in% c("exon_UTR", "lncRNA", "upstream_flank", "downstream_flank", "exon", "mRNA", "flanks", "nonCDS"),
          .(background, set, track, a, b, c, d, or = round(or_cmle, 3), lo = round(lo, 3), hi = round(hi, 3), p = signif(p, 3))])

cat("\n== shares (probability scale) ==\n")
for (s in c("ploidy DML", "heat DML")) for (t in c("gene", "intron", "CDS", "TE", "intergenic")) {
  x <- res[variant == "V1 all DML vs other background" & set == s & track == t]
  for (b in names(bgs)) {
    y <- x[background == b]
    cat(sprintf("%-10s %-10s %-18s DML in feature %.1f%%  other %.1f%%  risk ratio %.2f  | P(DML|in) %.3f%% P(DML|out) %.3f%% ratio %.2f\n",
                s, t, b, 100 * y$a / (y$a + y$b), 100 * y$c / (y$c + y$d), (y$a / (y$a + y$b)) / (y$c / (y$c + y$d)),
                100 * y$a / (y$a + y$c), 100 * y$b / (y$b + y$d), (y$a / (y$a + y$c)) / (y$b / (y$b + y$d))))
  }
}
saveRDS(cpg, file.path(scr, "cpg_features.rds"))
