# Table 1 loci recomputed from the 22 raw Bismark coverage files, without methylKit and without any pipeline output.
suppressPackageStartupMessages(library(data.table))
scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C7-legacy-go-semantics"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
samples <- readLines("/home/user/gannet_cache/cov_files/samples.txt")
meta <- data.table(sample = samples, heat = rep(c(1L, 0L), c(10, 12)), ploidy = c(rep(0L, 5), rep(1L, 5), rep(0L, 6), rep(1L, 6)))
# median coverage of the CpGs with at least 10 reads (what normalizeCoverage uses after filterByCoverage(lo.count = 10))
med <- sapply(samples, function(s) { h <- fread(file.path(scratch, "raw", paste0("hist_", s, ".txt")), col.names = c("cov", "n")); setorder(h, cov)
  n <- sum(h$n); cs <- cumsum(h$n); lo <- h$cov[which(cs >= floor((n + 1) / 2))[1]]; hi <- h$cov[which(cs >= ceiling((n + 1) / 2))[1]]; (lo + hi) / 2 })
meta[, median_cov := med[sample]][, sc := max(median_cov) / median_cov]
print(meta)
loci <- rbindlist(lapply(samples, function(s) { d <- fread(file.path(scratch, "raw", paste0("loci_", s, ".txt")), header = FALSE, fill = TRUE)
  d <- d[V2 != "LINES"]; setnames(d, c("sample", "line", "chr", "start", "end", "pct", "numC", "numT")); d[, line := as.integer(line)]; d }))
loci[, `:=`(start = as.integer(start), numC = as.numeric(numC), numT = as.numeric(numT))]
loci <- merge(loci, meta, by = "sample")
loci[, cov := numC + numT]
cat("loci x samples found:", nrow(loci), "of", 19 * 22, "; minimum coverage:", min(loci$cov), "\n")
# methylKit normalizeCoverage
loci[, ncov := round(sc * cov)][, `:=`(nC = round(ncov * numC / cov), nT = round(ncov * numT / cov))]
wm <- function(c, t) 100 * sum(c) / sum(c + t)
res <- loci[, .(
  heat_minus_control = wm(nC[heat == 1], nT[heat == 1]) - wm(nC[heat == 0], nT[heat == 0]),
  triploid_minus_diploid = wm(nC[ploidy == 1], nT[ploidy == 1]) - wm(nC[ploidy == 0], nT[ploidy == 0]),
  heat_effect_in_diploids = wm(nC[heat == 1 & ploidy == 0], nT[heat == 1 & ploidy == 0]) - wm(nC[heat == 0 & ploidy == 0], nT[heat == 0 & ploidy == 0]),
  heat_effect_in_triploids = wm(nC[heat == 1 & ploidy == 1], nT[heat == 1 & ploidy == 1]) - wm(nC[heat == 0 & ploidy == 1], nT[heat == 0 & ploidy == 1]),
  ploidy_effect_in_controls = wm(nC[heat == 0 & ploidy == 1], nT[heat == 0 & ploidy == 1]) - wm(nC[heat == 0 & ploidy == 0], nT[heat == 0 & ploidy == 0]),
  ploidy_effect_in_heat = wm(nC[heat == 1 & ploidy == 1], nT[heat == 1 & ploidy == 1]) - wm(nC[heat == 1 & ploidy == 0], nT[heat == 1 & ploidy == 0]),
  D_heat = wm(nC[heat == 1 & ploidy == 0], nT[heat == 1 & ploidy == 0]), T_heat = wm(nC[heat == 1 & ploidy == 1], nT[heat == 1 & ploidy == 1]),
  D_ctrl = wm(nC[heat == 0 & ploidy == 0], nT[heat == 0 & ploidy == 0]), T_ctrl = wm(nC[heat == 0 & ploidy == 1], nT[heat == 0 & ploidy == 1]),
  raw_heat_minus_control = wm(numC[heat == 1], numT[heat == 1]) - wm(numC[heat == 0], numT[heat == 0]),
  raw_triploid_minus_diploid = wm(numC[ploidy == 1], numT[ploidy == 1]) - wm(numC[ploidy == 0], numT[ploidy == 0])), by = .(chr, start)]
res[, difference_of_differences := heat_effect_in_triploids - heat_effect_in_diploids]
t1 <- fread(file.path(repo, "04_output/09_go_enrichment_legacy/tables/table_1_legacy.csv"))
m <- merge(t1, res, by.x = c("chromosome", "start"), by.y = c("chr", "start"))
m[, recomputed := fifelse(dml_set == "ploidy-DML", triploid_minus_diploid, heat_minus_control)]
m[, other_contrast := fifelse(dml_set == "ploidy-DML", heat_minus_control, triploid_minus_diploid)]
setorder(m, dml_set, gene_id, start)
options(width = 250)
cat("\nTable 1 value against the value recomputed from raw counts (normalized as methylKit does):\n")
print(m[, .(dml_set, gene = substr(gene_name, 1, 22), start, table1 = methylation_difference, recomputed = round(recomputed, 2),
            abs_diff = round(abs(methylation_difference - recomputed), 3), other_contrast = round(other_contrast, 2))])
cat("max |Table 1 - recomputed|:", max(abs(m$methylation_difference - m$recomputed)), "\n")
cat("\nThe three 'interaction' rows: group means (percent methylation, normalized counts pooled within group) and effects\n")
print(m[dml_set == "interaction-DML", .(gene = substr(gene_name, 1, 20), start, table1 = methylation_difference, D_ctrl = round(D_ctrl, 1), T_ctrl = round(T_ctrl, 1),
        D_heat = round(D_heat, 1), T_heat = round(T_heat, 1), heat_in_D = round(heat_effect_in_diploids, 1), heat_in_T = round(heat_effect_in_triploids, 1),
        diff_of_diff = round(difference_of_differences, 1), ploidy_in_ctrl = round(ploidy_effect_in_controls, 1), ploidy_in_heat = round(ploidy_effect_in_heat, 1))])
cat("\nPloidy rows: triploid minus diploid within each cohort\n")
print(m[dml_set == "ploidy-DML", .(gene = substr(gene_name, 1, 20), start, table1 = methylation_difference, all22 = round(triploid_minus_diploid, 1),
        in_controls = round(ploidy_effect_in_controls, 1), in_heat_cohort = round(ploidy_effect_in_heat, 1))])
cat("\nHeat rows: heat minus control within each ploidy\n")
print(m[dml_set == "heat-DML", .(gene = substr(gene_name, 1, 20), start, table1 = methylation_difference, all22 = round(heat_minus_control, 1),
        in_diploids = round(heat_effect_in_diploids, 1), in_triploids = round(heat_effect_in_triploids, 1))])
fwrite(m, file.path(scratch, "C_table1_from_raw.csv"))
# per-library values at the three interaction loci
cat("\nPer-library percent methylation (raw) at the three interaction loci\n")
x <- loci[paste(chr, start) %in% c("NC_047564.1 8425110", "NC_047561.1 28607016", "NW_022994979.1 42498")]
x[, grp := paste0(ifelse(ploidy == 1, "T", "D"), ifelse(heat == 1, "_heat", "_ctrl"))]
print(dcast(x[, .(v = paste0(round(100 * numC / cov), "% (", cov, ")")), by = .(chr, start, grp, sample)][, .(vals = paste(v, collapse = ", ")), by = .(chr, start, grp)], chr + start ~ grp, value.var = "vals"))
