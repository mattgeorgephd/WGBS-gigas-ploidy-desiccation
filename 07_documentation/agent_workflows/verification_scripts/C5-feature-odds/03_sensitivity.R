# Sensitivity of the odds ratios: eligible range, definition of the mean,
# stratification by methylation level and coverage, clustering, the 2.2 kb locus.
suppressPackageStartupMessages({library(data.table); library(splines); library(sandwich)})
scr  <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C5-feature-odds"
cpg <- readRDS(file.path(scr, "cpg_features.rds"))
N <- nrow(cpg)
set.seed(20261119)

or2 <- function(dml, bg, feat) {
  use <- bg | dml
  a <- sum(dml & feat); b <- sum(dml & !feat)
  cc <- sum(use & !dml & feat); d <- sum(use & !dml & !feat)
  ft <- fisher.test(matrix(c(a, cc, b, d), 2))
  list(a = a, b = b, c = cc, d = d, or = unname(ft$estimate), lo = ft$conf.int[1], hi = ft$conf.int[2], p = ft$p.value)
}
feats <- c("gene", "intron", "CDS", "TE", "intergenic", "exon_UTR", "lncRNA")
dsets <- c(ploidy = "dml_ploidy", heat = "dml_heat")

## ---- 0. probability scale ---------------------------------------------------
cat("== 0. shares on the probability scale (eligible = [10,90], mean of per-library %) ==\n")
elig <- cpg$mean_unw >= 10 & cpg$mean_unw <= 90
for (s in names(dsets)) for (f in c("gene", "intron", "CDS", "TE", "intergenic")) {
  dml <- cpg[[dsets[[s]]]]; ft <- cpg[[paste0("in_", f)]]
  for (bn in c("all", "eligible")) {
    bgv <- if (bn == "all") rep(TRUE, N) else elig
    x <- or2(dml, bgv, ft)
    cat(sprintf("%-6s %-10s vs %-8s: DML in feature %.1f%%, other CpGs %.1f%% (ratio of shares %.2f); P(DML|in) %.3f%%, P(DML|out) %.3f%% (ratio %.2f); OR %.2f\n",
                s, f, bn, 100 * x$a / (x$a + x$b), 100 * x$c / (x$c + x$d), (x$a / (x$a + x$b)) / (x$c / (x$c + x$d)),
                100 * x$a / (x$a + x$c), 100 * x$b / (x$b + x$d), (x$a / (x$a + x$c)) / (x$b / (x$b + x$d)), x$or))
  }
}

## ---- 1. eligible range --------------------------------------------------------
cat("\n== 1. odds ratios by eligible range ==\n")
ranges <- list(c(0, 100), c(2.5, 97.5), c(5, 95), c(10, 90), c(15, 85), c(20, 80), c(25, 75), c(30, 70))
out <- list()
for (def in c("mean_unw", "mean_pooled")) for (r in ranges) {
  m <- cpg[[def]]
  bgv <- m >= r[1] & m <= r[2]
  for (s in names(dsets)) for (f in feats) {
    dml <- cpg[[dsets[[s]]]]
    x <- or2(dml, bgv, cpg[[paste0("in_", f)]])
    x2 <- or2(dml & bgv, bgv, cpg[[paste0("in_", f)]])          # DML restricted to the range
    out[[length(out) + 1]] <- data.table(mean_def = def, lo_pct = r[1], hi_pct = r[2], n_bg = sum(bgv),
      set = s, n_dml = sum(dml), n_dml_in_range = sum(dml & bgv), track = f,
      or = x$or, ci_lo = x$lo, ci_hi = x$hi, p = x$p, or_dml_in_range_only = x2$or,
      lo_in_range_only = x2$lo, hi_in_range_only = x2$hi)
  }
}
rng <- rbindlist(out)
fwrite(rng, file.path(scr, "or_by_eligible_range.csv"))
show <- rng[mean_def == "mean_unw" & track %in% c("gene", "intron", "CDS", "TE", "intergenic")]
show[, lab := sprintf("%.2f (%.2f-%.2f)", or, ci_lo, ci_hi)]
print(dcast(show, set + track ~ paste0(lo_pct, "-", hi_pct), value.var = "lab")[
  , c("set", "track", "0-100", "2.5-97.5", "5-95", "10-90", "15-85", "20-80", "25-75", "30-70"), with = FALSE])
cat("\nbackground sizes and DML inside the range (mean of per-library %):\n")
print(unique(rng[mean_def == "mean_unw", .(lo_pct, hi_pct, n_bg, set, n_dml, n_dml_in_range)]))
cat("\nDML restricted to the range (drops DML outside it), mean of per-library %:\n")
show[, lab2 := sprintf("%.2f (%.2f-%.2f)", or_dml_in_range_only, lo_in_range_only, hi_in_range_only)]
print(dcast(show, set + track ~ paste0(lo_pct, "-", hi_pct), value.var = "lab2")[
  , c("set", "track", "5-95", "10-90", "15-85", "20-80", "25-75", "30-70"), with = FALSE])
cat("\nsame with the pooled mean (sum of C / sum of coverage):\n")
show <- rng[mean_def == "mean_pooled" & track %in% c("gene", "intron", "CDS", "TE", "intergenic")]
show[, lab := sprintf("%.2f (%.2f-%.2f)", or, ci_lo, ci_hi)]
print(dcast(show, set + track ~ paste0(lo_pct, "-", hi_pct), value.var = "lab")[
  , c("set", "track", "5-95", "10-90", "15-85"), with = FALSE])
cat("\nratio eligible OR / all-tested OR for gene body:\n")
g <- rng[mean_def == "mean_unw" & track == "gene"]
g[, or_all := or[lo_pct == 0], by = set]
print(g[, .(set, range = paste0(lo_pct, "-", hi_pct), or = round(or, 3), ratio_to_all = round(or / or_all, 3),
            log_ratio = round(log(or) / log(or_all), 3))])

## ---- 2. heterogeneity inside the eligible band ---------------------------------
cat("\n== 2. P(DML) and feature composition by mean methylation bin (all tested CpGs) ==\n")
cpg[, mbin := cut(mean_unw, c(0, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 85, 90, 95, 100), include.lowest = TRUE)]
print(cpg[, .(n = .N, ploidy_dml = sum(dml_ploidy), heat_dml = sum(dml_heat),
              p_ploidy_pct = round(100 * mean(dml_ploidy), 3), p_heat_pct = round(100 * mean(dml_heat), 3),
              pct_gene = round(100 * mean(in_gene), 1), pct_intron = round(100 * mean(in_intron), 1),
              pct_CDS = round(100 * mean(in_CDS), 1), pct_TE = round(100 * mean(in_TE), 1),
              pct_intergenic = round(100 * mean(in_intergenic), 1), mean_cov = round(mean(mean_cov), 1)), keyby = mbin])

cat("\nmean coverage and mean methylation of eligible CpGs by feature:\n")
for (f in c("gene", "intron", "CDS", "TE", "intergenic", "exon_UTR", "lncRNA")) {
  i <- elig & cpg[[paste0("in_", f)]]; j <- elig & !cpg[[paste0("in_", f)]]
  cat(sprintf("%-11s in: n %6d  median meth %.1f  median mean-cov %.1f | out: n %6d  median meth %.1f  median mean-cov %.1f\n",
              f, sum(i), median(cpg$mean_unw[i]), median(cpg$mean_cov[i]), sum(j), median(cpg$mean_unw[j]), median(cpg$mean_cov[j])))
}

## ---- 3. Mantel-Haenszel within the eligible band ---------------------------------
cat("\n== 3. Mantel-Haenszel odds ratios among eligible CpGs, stratified ==\n")
e <- cpg[elig | dml_ploidy | dml_heat]
e[, m8 := cut(mean_unw, seq(0, 100, 10), include.lowest = TRUE)]
e[, m5 := cut(mean_unw, seq(0, 100, 5), include.lowest = TRUE)]
e[, covq := cut(mean_cov, quantile(mean_cov, 0:4 / 4), include.lowest = TRUE)]
mh <- function(d, y, f, strata) {
  tb <- table(d[[y]], d[[f]], interaction(d[, strata, with = FALSE], drop = TRUE))
  keep <- apply(tb, 3, function(z) all(rowSums(z) > 0) && all(colSums(z) > 0))
  r <- mantelhaen.test(tb[, , keep, drop = FALSE])
  c(or = unname(r$estimate), lo = r$conf.int[1], hi = r$conf.int[2], p = r$p.value, n_strata = sum(keep))
}
res3 <- list()
for (s in names(dsets)) for (f in feats) {
  y <- dsets[[s]]; fc <- paste0("in_", f)
  d <- e[elig | get(y)]
  crude <- or2(d[[y]], rep(TRUE, nrow(d)), d[[fc]])
  r1 <- mh(d, y, fc, "m8"); r2 <- mh(d, y, fc, "m5"); r3 <- mh(d, y, fc, c("m5", "covq"))
  res3[[length(res3) + 1]] <- data.table(set = s, track = f, crude = crude$or,
    mh_meth10 = r1["or"], mh_meth10_lo = r1["lo"], mh_meth10_hi = r1["hi"],
    mh_meth5 = r2["or"], mh_meth5_lo = r2["lo"], mh_meth5_hi = r2["hi"],
    mh_meth5_cov4 = r3["or"], mh_meth5_cov4_lo = r3["lo"], mh_meth5_cov4_hi = r3["hi"], p_meth5_cov4 = r3["p"])
}
res3 <- rbindlist(res3)
fwrite(res3, file.path(scr, "or_mantel_haenszel.csv"))
print(res3[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 3) else x)])

## ---- 4. logistic regression with cluster-robust standard errors -------------------
cat("\n== 4. logistic regression among eligible CpGs: feature + ns(mean methylation, 6) + log(mean coverage) ==\n")
e[, win10k := paste(chr, start %/% 10000L)]
e[, win100k := paste(chr, start %/% 100000L)]
res4 <- list()
for (s in names(dsets)) for (f in feats) {
  y <- dsets[[s]]; fc <- paste0("in_", f)
  d <- e[elig | get(y)]
  d[, yy := as.integer(get(y))][, ff := as.integer(get(fc))]
  m0 <- glm(yy ~ ff, binomial, d)
  m1 <- glm(yy ~ ff + ns(mean_unw, 6) + log(mean_cov), binomial, d)
  se <- function(m, cl) sqrt(diag(vcovCL(m, cluster = d[[cl]], type = "HC0"))["ff"])
  b0 <- coef(m0)["ff"]; b1 <- coef(m1)["ff"]
  s0n <- sqrt(diag(vcov(m0))["ff"]); s0a <- se(m0, "win10k"); s0b <- se(m0, "win100k")
  s1n <- sqrt(diag(vcov(m1))["ff"]); s1a <- se(m1, "win10k"); s1b <- se(m1, "win100k")
  res4[[length(res4) + 1]] <- data.table(set = s, track = f,
    or_crude = exp(b0), crude_naive = sprintf("%.2f-%.2f", exp(b0 - 1.96 * s0n), exp(b0 + 1.96 * s0n)),
    crude_cl10k = sprintf("%.2f-%.2f", exp(b0 - 1.96 * s0a), exp(b0 + 1.96 * s0a)),
    crude_cl100k = sprintf("%.2f-%.2f", exp(b0 - 1.96 * s0b), exp(b0 + 1.96 * s0b)),
    p_crude_cl100k = 2 * pnorm(-abs(b0 / s0b)),
    or_adj = exp(b1), adj_naive = sprintf("%.2f-%.2f", exp(b1 - 1.96 * s1n), exp(b1 + 1.96 * s1n)),
    adj_cl100k = sprintf("%.2f-%.2f", exp(b1 - 1.96 * s1b), exp(b1 + 1.96 * s1b)),
    p_adj_cl100k = 2 * pnorm(-abs(b1 / s1b)), se_inflation_100k = s0b / s0n)
}
res4 <- rbindlist(res4)
fwrite(res4, file.path(scr, "or_logistic_cluster.csv"))
print(res4[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 3) else x)])

## ---- 5. block bootstrap of the crude 2x2 odds ratios --------------------------------
cat("\n== 5. block bootstrap (blocks = 100 kb windows; 2000 resamples) of the crude odds ratios ==\n")
boot_or <- function(d, y, fc, block, B = 2000) {
  x <- d[, .(a = sum(get(y) & get(fc)), b = sum(get(y) & !get(fc)),
             c = sum(!get(y) & get(fc)), d = sum(!get(y) & !get(fc))), by = block]
  M <- as.matrix(x[, .(a, b, c, d)]); nb <- nrow(M)
  ors <- numeric(B)
  for (i in seq_len(B)) {
    w <- tabulate(sample.int(nb, nb, replace = TRUE), nb)
    t <- colSums(M * w)
    ors[i] <- (t[1] * t[4]) / (t[2] * t[3])
  }
  c(nblocks = nb, quantile(ors, c(0.025, 0.975), na.rm = TRUE))
}
res5 <- list()
for (bgname in c("all", "eligible")) for (s in names(dsets)) for (f in c("gene", "intron", "CDS", "TE", "intergenic")) {
  y <- dsets[[s]]; fc <- paste0("in_", f)
  d <- if (bgname == "all") copy(cpg) else cpg[elig | get(y)]
  d[, win100k := paste(chr, start %/% 100000L)]
  cr <- or2(d[[y]], rep(TRUE, nrow(d)), d[[fc]])
  b <- boot_or(d, y, fc, "win100k")
  res5[[length(res5) + 1]] <- data.table(background = bgname, set = s, track = f, or = cr$or,
    fisher_ci = sprintf("%.2f-%.2f", cr$lo, cr$hi), nblocks = b[1], boot_ci = sprintf("%.2f-%.2f", b[2], b[3]),
    boot_lo = b[2], boot_hi = b[3])
}
res5 <- rbindlist(res5)
fwrite(res5, file.path(scr, "or_block_bootstrap.csv"))
print(res5[, .(background, set, track, or = round(or, 3), fisher_ci, boot_ci, nblocks)])

## ---- 6. the 2.2 kb locus on NC_047563.1 ------------------------------------------------
cat("\n== 6. without NC_047563.1:39,633,323-39,635,520 ==\n")
loc <- cpg$chr == "NC_047563.1" & cpg$start >= 39633323 - 1 & cpg$start <= 39635520
cat("tested CpGs in locus:", sum(loc), " eligible:", sum(loc & elig), " ploidy DML:", sum(loc & cpg$dml_ploidy),
    " heat DML:", sum(loc & cpg$dml_heat), " shared:", sum(loc & cpg$dml_shared), "\n")
cat("locus CpGs by track: ", paste(sapply(c("gene", "CDS", "exon_UTR", "intron", "lncRNA", "TE", "intergenic"),
    function(f) sprintf("%s %d", f, sum(loc & cpg[[paste0("in_", f)]]))), collapse = "; "), "\n")
for (s in names(dsets)) for (f in feats) {
  dml <- cpg[[dsets[[s]]]] & !loc
  x <- or2(dml, elig & !loc, cpg[[paste0("in_", f)]]); x0 <- or2(cpg[[dsets[[s]]]], elig, cpg[[paste0("in_", f)]])
  cat(sprintf("%-6s %-10s eligible: with locus %.3f (%.2f-%.2f, p=%.2g) ; without locus %.3f (%.2f-%.2f, p=%.2g)\n",
              s, f, x0$or, x0$lo, x0$hi, x0$p, x$or, x$lo, x$hi, x$p))
}

## ---- 7. within gene bodies: CDS against intron ------------------------------------------
cat("\n== 7. among eligible CpGs inside genes: CDS-only vs intron-only ==\n")
for (s in names(dsets)) {
  dml <- cpg[[dsets[[s]]]]
  g <- (elig | dml) & cpg$in_gene & xor(cpg$in_CDS, cpg$in_intron)
  tb <- table(dml = dml[g], CDS = cpg$in_CDS[g])
  ft <- fisher.test(tb)
  cat(sprintf("%-6s n genic eligible in exactly one of CDS/intron %d ; P(DML|CDS) %.3f%% P(DML|intron) %.3f%% ; OR %.2f (%.2f-%.2f)\n",
              s, sum(g), 100 * tb["TRUE", "TRUE"] / sum(tb[, "TRUE"]), 100 * tb["TRUE", "FALSE"] / sum(tb[, "FALSE"]),
              ft$estimate, ft$conf.int[1], ft$conf.int[2]))
}
cat("\ntrack co-occurrence among all tested CpGs (n):\n")
tt <- c("gene", "CDS", "exon_UTR", "intron", "lncRNA", "TE", "intergenic", "upstream_flank", "downstream_flank")
co <- sapply(tt, function(a) sapply(tt, function(b) sum(cpg[[paste0("in_", a)]] & cpg[[paste0("in_", b)]])))
print(co)
cat("CpGs in none of the eight legacy tracks:", sum(!Reduce(`|`, lapply(c("CDS", "exon_UTR", "intron", "lncRNA", "upstream_flank", "downstream_flank", "intergenic", "TE"), function(f) cpg[[paste0("in_", f)]]))), "\n")
cat("CpGs in neither gene nor intergenic nor flanks:", sum(!cpg$in_gene & !cpg$in_intergenic & !cpg$in_flanks), "\n")
