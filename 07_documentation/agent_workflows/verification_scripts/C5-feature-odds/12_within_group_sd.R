# Does DML location among eligible CpGs track between-animal variability within groups?
suppressPackageStartupMessages({library(data.table); library(splines); library(sandwich)})
scr <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C5-feature-odds"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
cpg <- readRDS(file.path(scr, "cpg_features.rds"))
u <- fread(file.path(repo, "04_output/05_differential_methylation/tables/united_cpg_legacy.csv.gz"))
stopifnot(identical(u$chr, cpg$chr), identical(u$start, cpg$start))
nc <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]); nt <- as.matrix(u[, paste0("numTs", 1:22), with = FALSE])
pm <- 100 * nc / (nc + nt); rm(u, nc, nt); invisible(gc())
grp <- list(dip_heat = 1:5, tri_heat = 6:10, dip_ctrl = 11:16, tri_ctrl = 17:22)   # legacy_sample_index order
rv <- function(m) { k <- ncol(m); (rowSums(m^2) - k * rowMeans(m)^2) / (k - 1) }
v <- sapply(grp, function(i) rv(pm[, i, drop = FALSE]))
df <- sapply(grp, length) - 1
cpg[, sd_within := sqrt(as.vector(v %*% df) / sum(df))]         # pooled within-group SD, percentage points
rm(pm, v); invisible(gc())
e <- cpg[elig_p == TRUE | dml_ploidy | dml_heat]
e[, m10 := cut(mean_ct, seq(0, 100, 10), include.lowest = TRUE)]
cat("== median pooled within-group SD (percentage points) of eligible CpGs, by feature and methylation bin ==\n")
print(e[elig_p == TRUE, .(n = .N, all = round(median(sd_within), 1), gene = round(median(sd_within[in_gene]), 1), nongene = round(median(sd_within[!in_gene]), 1),
        CDS = round(median(sd_within[in_CDS]), 1), intron = round(median(sd_within[in_intron]), 1), TE = round(median(sd_within[in_TE]), 1),
        intergenic = round(median(sd_within[in_intergenic]), 1)), keyby = m10])
cat("\n== P(DML) by quartile of within-group SD among eligible CpGs ==\n")
e[, sdq := cut(sd_within, quantile(sd_within, 0:4 / 4), include.lowest = TRUE)]
print(e[, .(n = .N, p_ploidy_pct = round(100 * mean(dml_ploidy), 2), p_heat_pct = round(100 * mean(dml_heat), 2), pct_gene = round(100 * mean(in_gene), 1),
            pct_CDS = round(100 * mean(in_CDS), 1), pct_TE = round(100 * mean(in_TE), 1), median_cov = round(median(mean_cov), 1)), keyby = sdq])
e[, win := paste(chr, start %/% 100000L)]
cat("\n== feature odds ratio among eligible CpGs, adding the within-group SD (clustered 95% intervals, 100 kb) ==\n")
for (s in c("dml_ploidy", "dml_heat")) for (f in c("gene", "intron", "CDS", "TE", "intergenic")) {
  d <- e[elig_p == TRUE | get(s)]; d[, yy := as.integer(get(s))][, ff := as.integer(get(paste0("in_", f)))]
  out <- sapply(list(yy ~ ff, yy ~ ff + ns(mean_ct, 6), yy ~ ff + ns(mean_ct, 6) + ns(sd_within, 4)), function(fm) {
    m <- glm(fm, binomial, d); b <- coef(m)["ff"]; sc <- sqrt(diag(vcovCL(m, cluster = d$win, type = "HC0"))["ff"])
    sprintf("%.2f (%.2f-%.2f)", exp(b), exp(b - 1.96 * sc), exp(b + 1.96 * sc)) })
  cat(sprintf("%-10s %-10s crude %s | +meth %s | +meth +within-group SD %s\n", s, f, out[1], out[2], out[3]))
}
