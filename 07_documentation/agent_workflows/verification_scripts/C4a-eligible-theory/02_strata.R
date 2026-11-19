# Check (b)/(c): within the eligible set, how uneven is the chance of being a DML,
# and do the location odds ratios survive stratification on mean methylation?
suppressPackageStartupMessages({library(data.table); library(GenomicRanges); library(splines)})
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
sc <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4a-eligible-theory"
mode <- "legacy"
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv"))
trip <- lib$ploidy == "triploid"; heat <- lib$treatment == "heat"
u <- fread(file.path(repo, sprintf("04_output/05_differential_methylation/tables/united_cpg_%s.csv.gz", mode)))
C <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]); Tt <- as.matrix(u[, paste0("numTs", 1:22), with = FALSE])
cov <- C + Tt; pm <- 100 * C / cov
cpg <- u[, .(chr, start)]; rm(u, C, Tt); gc(verbose = FALSE)
m <- rowMeans(pm)
cell <- interaction(lib$ploidy, lib$treatment, drop = TRUE)
# within-cell residual SD (pooled over the four ploidy-by-cohort cells): variability not due to either contrast
res <- pm; for (l in levels(cell)) res[, cell == l] <- pm[, cell == l] - rowMeans(pm[, cell == l, drop = FALSE])
sd_within <- sqrt(rowSums(res^2) / (22 - 4)); rm(res)
dml <- list(ploidy = fread(file.path(repo, sprintf("04_output/05_differential_methylation/tables/dml_ploidy_%s.csv", mode))),
            heat   = fread(file.path(repo, sprintf("04_output/05_differential_methylation/tables/dml_heat_%s.csv", mode))))
cpg[, `:=`(mean_pm = m, sd_within = sd_within, mean_cov = rowMeans(cov), dml_ploidy = seq_len(.N) %in% dml$ploidy$united_row, dml_heat = seq_len(.N) %in% dml$heat$united_row,
           m_heat = rowMeans(pm[, heat]), m_ctrl = rowMeans(pm[, !heat]), m_trip = rowMeans(pm[, trip]), m_dip = rowMeans(pm[, !trip]))]
cpg[, eligible := mean_pm >= 10 & mean_pm <= 90]
gr <- GRanges(cpg$chr, IRanges(cpg$start + 1L, cpg$start + 2L))
rd <- function(file, bed = FALSE) {
  p <- file.path(repo, "03_input_files/genome_features", paste0("cgigas_uk_roslin_v1_", file, ".gz"))
  if (!bed) d <- fread(p, header = FALSE, sep = "\t", quote = "", select = c(1, 4, 5), col.names = c("chr", "start", "end"))
  else { d <- fread(p, header = FALSE, sep = "\t", quote = "", select = 1:3, col.names = c("chr", "start", "end")); d[, start := start + 1L] }
  d <- d[end >= start]; GRanges(d$chr, IRanges(d$start, d$end))
}
tracks <- list(gene = rd("gene.gff"), CDS = rd("CDS.gff"), intron = rd("intron.bed", TRUE), TE = rd("rm.te.bed", TRUE), intergenic = rd("intergenic.bed", TRUE))
for (t in names(tracks)) cpg[, (t) := suppressWarnings(overlapsAny(gr, tracks[[t]], ignore.strand = TRUE))]
saveRDS(cpg, file.path(sc, "cpg_table_legacy.rds"))

el <- cpg[eligible == TRUE | dml_ploidy | dml_heat]
cat(sprintf("eligible (plus the 2 outside DML): %d\n", nrow(el)))
el[, mean_bin := cut(mean_pm, c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100), right = FALSE)]
cat("\nP(DML) by mean-methylation bin within the eligible set:\n")
print(el[, .(n = .N, ploidy_dml = sum(dml_ploidy), pct_ploidy = round(100 * mean(dml_ploidy), 2), heat_dml = sum(dml_heat), pct_heat = round(100 * mean(dml_heat), 2),
             pct_in_gene = round(100 * mean(gene), 1), pct_CDS = round(100 * mean(CDS), 1), pct_intron = round(100*mean(intron), 1), pct_TE = round(100 * mean(TE), 1), pct_intergenic = round(100 * mean(intergenic), 1),
             median_sd_within = round(median(sd_within), 1)), by = mean_bin][order(mean_bin)])
el[, sd_bin := cut(sd_within, c(0, 5, 10, 15, 20, 25, 30, 100), right = FALSE)]
cat("\nP(DML) by within-cell SD bin (eligible set):\n")
print(el[, .(n = .N, pct_ploidy = round(100 * mean(dml_ploidy), 2), pct_heat = round(100 * mean(dml_heat), 2), pct_in_gene = round(100 * mean(gene), 1), pct_CDS = round(100*mean(CDS), 1), pct_TE = round(100*mean(TE),1), pct_intergenic = round(100 * mean(intergenic), 1)), by = sd_bin][order(sd_bin)])

mh <- function(y, x, strata) {
  tb <- table(factor(y, c(TRUE, FALSE)), factor(x, c(TRUE, FALSE)), droplevels(factor(strata)))
  keep <- apply(tb, 3, function(z) all(rowSums(z) > 0) && all(colSums(z) > 0))
  r <- mantelhaen.test(tb[, , keep, drop = FALSE], correct = FALSE)
  c(or = unname(r$estimate), lo = r$conf.int[1], hi = r$conf.int[2], p = r$p.value)
}
crude <- function(y, x) { ft <- fisher.test(table(factor(y, c(TRUE, FALSE)), factor(x, c(TRUE, FALSE)))); c(or = unname(ft$estimate), lo = ft$conf.int[1], hi = ft$conf.int[2], p = ft$p.value) }
el[, mean_bin5 := cut(mean_pm, seq(0, 100, 5), right = FALSE)]
el[, sd_bin5 := cut(sd_within, c(0, 4, 8, 12, 16, 20, 25, 100), right = FALSE)]
el[, cov_bin := cut(mean_cov, quantile(mean_cov, 0:5/5), include.lowest = TRUE)]
out <- list()
for (s in c("ploidy", "heat")) for (t in names(tracks)) {
  y <- el[[paste0("dml_", s)]]; x <- el[[t]]
  out[[length(out) + 1]] <- data.table(set = s, track = t, analysis = c("crude (eligible background)", "MH, strata = mean (5-pt bins)", "MH, strata = mean x within-cell SD", "MH, strata = mean x SD x coverage quintile"),
     rbind(crude(y, x), mh(y, x, el$mean_bin5), mh(y, x, interaction(el$mean_bin5, el$sd_bin5, drop = TRUE)), mh(y, x, interaction(el$mean_bin5, el$sd_bin5, el$cov_bin, drop = TRUE))))
}
out <- rbindlist(out); out[, `:=`(or = round(or, 2), lo = round(lo, 2), hi = round(hi, 2), p = signif(p, 3))]
cat("\nLocation odds ratios within the eligible set, crude and stratified:\n"); print(out)
fwrite(out, file.path(sc, "location_or_stratified.csv"))

# logistic regression with smooth terms
for (s in c("ploidy", "heat")) {
  f0 <- glm(as.formula(sprintf("dml_%s ~ gene", s)), binomial, el)
  f1 <- glm(as.formula(sprintf("dml_%s ~ gene + ns(mean_pm, 6)", s)), binomial, el)
  f2 <- glm(as.formula(sprintf("dml_%s ~ gene + ns(mean_pm, 6) + ns(sd_within, 5) + ns(log(mean_cov), 4)", s)), binomial, el)
  cat(sprintf("[%s] logistic OR for gene body: crude %.2f; + spline(mean) %.2f; + spline(mean) + spline(within-cell SD) + spline(log coverage) %.2f\n", s, exp(coef(f0)["geneTRUE"]), exp(coef(f1)["geneTRUE"]), exp(coef(f2)["geneTRUE"])))
  for (t in c("CDS", "intron", "TE", "intergenic")) {
    g0 <- glm(as.formula(sprintf("dml_%s ~ %s", s, t)), binomial, el)
    g2 <- glm(as.formula(sprintf("dml_%s ~ %s + ns(mean_pm, 6) + ns(sd_within, 5) + ns(log(mean_cov), 4)", s, t)), binomial, el)
    g1 <- glm(as.formula(sprintf("dml_%s ~ %s + ns(mean_pm, 6)", s, t)), binomial, el)
    cat(sprintf("    %-10s crude %.2f; + mean %.2f; + mean + SD + coverage %.2f\n", t, exp(coef(g0)[2]), exp(coef(g1)[2]), exp(coef(g2)[2])))
  }
}

# Counterfactual eligibility: is a DML eligible only because of the difference itself?
cat("\nGroup means of DML (unweighted):\n")
h <- cpg[dml_heat == TRUE]; p <- cpg[dml_ploidy == TRUE]
cat(sprintf("heat DML %d: control-group mean <10: %d; >90: %d; heat-group mean <10: %d; >90: %d; lower of the two group means <10: %d; higher of the two >90: %d; either: %d\n",
            nrow(h), sum(h$m_ctrl < 10), sum(h$m_ctrl > 90), sum(h$m_heat < 10), sum(h$m_heat > 90), sum(pmin(h$m_ctrl, h$m_heat) < 10), sum(pmax(h$m_ctrl, h$m_heat) > 90),
            sum(pmin(h$m_ctrl, h$m_heat) < 10 | pmax(h$m_ctrl, h$m_heat) > 90)))
cat(sprintf("ploidy DML %d: diploid mean <10: %d; >90: %d; triploid mean <10: %d; >90: %d; lower of two <10: %d; higher of two >90: %d; either: %d\n",
            nrow(p), sum(p$m_dip < 10), sum(p$m_dip > 90), sum(p$m_trip < 10), sum(p$m_trip > 90), sum(pmin(p$m_dip, p$m_trip) < 10), sum(pmax(p$m_dip, p$m_trip) > 90),
            sum(pmin(p$m_dip, p$m_trip) < 10 | pmax(p$m_dip, p$m_trip) > 90)))
# Background comparison: non-DML eligible CpGs
nb <- cpg[eligible == TRUE & !dml_heat & !dml_ploidy]
cat(sprintf("non-DML eligible %d: lower of heat/control group means <10: %d (%.1f%%); higher >90: %d (%.1f%%)\n", nrow(nb), sum(pmin(nb$m_ctrl, nb$m_heat) < 10), 100*mean(pmin(nb$m_ctrl, nb$m_heat) < 10),
            sum(pmax(nb$m_ctrl, nb$m_heat) > 90), 100 * mean(pmax(nb$m_ctrl, nb$m_heat) > 90)))
