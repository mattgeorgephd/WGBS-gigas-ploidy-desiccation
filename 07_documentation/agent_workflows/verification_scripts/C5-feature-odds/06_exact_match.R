# (i) eligible set with the pipeline's definition of percent methylation, numCs/(numCs+numTs);
# (ii) exact comparison of my tables with the pipeline tables; (iii) workbook counts; (iv) linkage groups.
suppressPackageStartupMessages({library(data.table); library(readxl)})
scr  <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C5-feature-odds"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
cpg <- readRDS(file.path(scr, "cpg_features.rds"))
u <- fread(file.path(repo, "04_output/05_differential_methylation/tables/united_cpg_legacy.csv.gz"))
nc <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]); nt <- as.matrix(u[, paste0("numTs", 1:22), with = FALSE])
cv <- as.matrix(u[, paste0("coverage", 1:22), with = FALSE])
cat("cells where coverage != numCs + numTs:", sum(cv != nc + nt), "of", length(cv),
    sprintf("(%.2f%%); range of difference: %s\n", 100 * mean(cv != nc + nt), paste(range(cv - nc - nt), collapse = " to ")))
cpg[, mean_ct := rowMeans(100 * nc / (nc + nt))]
rm(u, nc, nt, cv); invisible(gc())
cat("eligible, numCs/(numCs+numTs): ", sum(cpg$mean_ct >= 10 & cpg$mean_ct <= 90), "\n")
cat("eligible, numCs/coverage:      ", sum(cpg$mean_unw >= 10 & cpg$mean_unw <= 90), "\n")
cat("in one set but not the other:  ", sum(xor(cpg$mean_ct >= 10 & cpg$mean_ct <= 90, cpg$mean_unw >= 10 & cpg$mean_unw <= 90)), "\n")
cat("below 10 / above 90 (pipeline definition): ", sum(cpg$mean_ct < 10), sprintf("(%.2f%%)", 100 * mean(cpg$mean_ct < 10)), "/",
    sum(cpg$mean_ct > 90), sprintf("(%.2f%%)", 100 * mean(cpg$mean_ct > 90)), "\n")
cpg[, elig_p := mean_ct >= 10 & mean_ct <= 90]
cat("DML outside (pipeline definition): ploidy", sum(cpg$dml_ploidy & !cpg$elig_p), " heat", sum(cpg$dml_heat & !cpg$elig_p), "\n")
N <- nrow(cpg)

sets <- list("all analysed CpGs" = rep(TRUE, N), "eligible CpGs" = cpg$elig_p, "ploidy DML" = cpg$dml_ploidy,
  "heat DML" = cpg$dml_heat, "shared" = cpg$dml_shared, "shared, both hyper" = cpg$shared_class %in% "both hyper",
  "shared, both hypo" = cpg$shared_class %in% "both hypo", "shared, opposite" = cpg$shared_class %in% "opposite")
tr <- sub("^in_", "", grep("^in_", names(cpg), value = TRUE))
eight <- c("CDS", "exon_UTR", "intron", "lncRNA", "upstream_flank", "downstream_flank", "intergenic", "TE")
ov <- rbindlist(lapply(names(sets), function(s) {
  i <- sets[[s]]
  d <- data.table(set = s, n_loci = sum(i), track = tr, n_overlap = vapply(tr, function(t) sum(cpg[[paste0("in_", t)]][i]), numeric(1)))
  d[, pct_of_loci := 100 * n_overlap / n_loci]
  d[track %in% eight, share8 := 100 * n_overlap / sum(n_overlap)]
  d
}))
pipe <- fread(file.path(repo, "04_output/07_genome_features/tables/feature_overlap_by_track_legacy.csv"))
cmp <- merge(ov, pipe, by = c("set", "track"), suffixes = c("_mine", "_pipe"))
cat("\n== overlap table, pipeline eligible definition: rows", nrow(cmp), "\n")
cat("n_loci equal:", all(cmp$n_loci_mine == cmp$n_loci_pipe), " n_overlap equal:", all(cmp$n_overlap_mine == cmp$n_overlap_pipe),
    " max|share diff|:", max(abs(cmp$share8 - cmp$legacy_share_pct), na.rm = TRUE),
    " max|pct_of_loci diff|:", max(abs(cmp$pct_of_loci_mine - cmp$pct_of_loci_pipe)), "\n")
fwrite(ov, file.path(scr, "my_overlap_by_track_pipeline_eligible.csv"))

## Fisher tests exactly as the claim's table
pa <- fread(file.path(repo, "04_output/07_genome_features/tables/feature_association_tests_legacy.csv"))
mine <- rbindlist(lapply(c("all analysed CpGs", "eligible CpGs"), function(b) rbindlist(lapply(c("ploidy DML", "heat DML", "shared"), function(s)
  rbindlist(lapply(c(eight, "gene"), function(t) {
    dml <- sets[[s]]; bgv <- sets[[b]]; f <- cpg[[paste0("in_", t)]]
    a <- sum(dml & f); b2 <- sum(dml & !f); cc <- sum(bgv & !dml & f); d <- sum(bgv & !dml & !f)
    ft <- fisher.test(matrix(c(a, cc, b2, d), 2))
    data.table(background = b, set = s, track = t, a = a, dml_n = a + b2, other_in = cc, other_n = cc + d,
               or = unname(ft$estimate), lo = ft$conf.int[1], hi = ft$conf.int[2], p = ft$p.value)
  }))))))
mine[, q := p.adjust(p, "BH")]
m <- merge(mine, pa, by = c("background", "set", "track"))
cat("\n== association table: rows mine", nrow(mine), " pipeline", nrow(pa), " merged", nrow(m), "\n")
cat("dml_in_track equal:", all(m$a == m$dml_in_track), " dml_total equal:", all(m$dml_n == m$dml_total),
    "\nmax |pct_other diff|:", max(abs(100 * m$other_in / m$other_n - m$pct_other_cpg)),
    "\nmax rel OR diff:", max(abs(m$or / m$odds_ratio - 1), na.rm = TRUE),
    "\nmax rel CI diff:", max(abs(m$lo / m$ci_low - 1), abs(m$hi / m$ci_high - 1), na.rm = TRUE),
    "\nmax rel p diff:", max(abs(m$p / m$p_fisher - 1)), "\nmax rel q diff:", max(abs(m$q / m$q_BH - 1)), "\n")
cat("OR of 0 (both):", sum(m$or == 0 & m$odds_ratio == 0), "\n")

## legacy workbook: DML overlap counts
f <- file.path(repo, "08_archive/analyses/DML_locations-all-datasets.xlsx")
d <- suppressMessages(as.data.table(read_excel(f, sheet = "data", col_names = FALSE, .name_repair = "minimal")))
setnames(d, paste0("V", seq_len(ncol(d))))
wb <- d[c(3:10, 12:16), .(feature = V1, bg10x = as.numeric(V10), ploidy = as.numeric(V11), heat = as.numeric(V12), shared = as.numeric(V13),
                         hypo = as.numeric(V14), hyper = as.numeric(V15), contrast = as.numeric(V16))]
wb[, track := c("exon_UTR", "CDS", "intron", "upstream_flank", "downstream_flank", "lncRNA", "TE", "intergenic", "gene", "exon", "mRNA", "nonCDS", "flanks")]
map <- c(ploidy = "ploidy DML", heat = "heat DML", shared = "shared", hypo = "shared, both hypo", hyper = "shared, both hyper", contrast = "shared, opposite")
wl <- melt(wb, id.vars = c("feature", "track"), measure.vars = names(map), variable.name = "col", value.name = "n_workbook")
wl[, set := map[as.character(col)]]
w <- merge(wl, ov, by = c("set", "track"))
cat("\n== legacy workbook DML counts vs my overlap counts: cells", nrow(w), " equal:", sum(w$n_workbook == w$n_overlap), "\n")
print(w[n_workbook != n_overlap])
bg <- merge(wb[, .(track, bg10x)], ov[set == "all analysed CpGs", .(track, n_overlap)], by = "track")
bg[, diff := bg10x - n_overlap][, ratio := bg10x / n_overlap]
cat("\nlegacy 10x background counts vs the tested CpGs:\n"); print(bg)
b8 <- bg[track %in% eight]
b8[, legacy_share := 100 * bg10x / sum(bg10x)][, regenerated_share := 100 * n_overlap / sum(n_overlap)][, d := regenerated_share - legacy_share]
cat("\nshares of the eight tracks, legacy background vs tested CpGs (sum legacy", sum(b8$bg10x), ", sum tested", sum(b8$n_overlap), "):\n")
print(b8[, .(track, legacy_share = round(legacy_share, 3), regenerated_share = round(regenerated_share, 3), diff = round(d, 3))])
cat("max |difference|:", max(abs(b8$d)), "\n")

## manuscript percentages
ms <- data.table(
  set = c(rep("all analysed CpGs", 8), "ploidy DML", "heat DML", "ploidy DML", "heat DML", "ploidy DML", "heat DML", "ploidy DML", "heat DML",
          "shared", "shared", "shared", "shared", "shared, both hypo", "shared, both hyper", rep("shared, opposite", 4)),
  track = c("intron", "CDS", "intergenic", "TE", "exon_UTR", "upstream_flank", "lncRNA", "downstream_flank",
            "intron", "intron", "CDS", "CDS", "intergenic", "intergenic", "TE", "TE",
            "intergenic", "TE", "lncRNA", "exon_UTR", "intron", "intron", "lncRNA", "exon_UTR", "TE", "intergenic"),
  manuscript = c(27.7, 21.0, 18.5, 18.4, 7.0, 2.8, 2.5, 2.1, 39.6, 42.0, 28.6, 27.6, 2.6, 3.2, 14.9, 16.0,
                 1.5, 10.3, 15.5, 11.1, 47.5, 43.2, 23.4, 14.9, 5.7, 1.1))
ms <- merge(ms, ov[, .(set, track, mine = share8)], by = c("set", "track"), sort = FALSE)
ms <- merge(ms, pipe[, .(set, track, pipeline = legacy_share_pct)], by = c("set", "track"), sort = FALSE)
ms[, match_1dp := abs(round(pipeline, 1) - manuscript) < 1e-9][, diff := round(pipeline - manuscript, 3)]
cat("\n== manuscript Results 3.3 percentages vs legacy_share_pct ==\n")
print(ms[, .(set, track, manuscript, pipeline = round(pipeline, 4), mine = round(mine, 4), diff, match_1dp)])
cat("DML percentages matching at one decimal:", sum(ms[set != "all analysed CpGs"]$match_1dp), "of", nrow(ms[set != "all analysed CpGs"]),
    "; background percentages matching:", sum(ms[set == "all analysed CpGs"]$match_1dp), "of 8\n")

## linkage groups
lg <- sprintf("NC_0475%02d.1", 59:68)
cpg[, lgi := match(chr, lg)]
tab <- cpg[, .(tested = .N, ploidy = sum(dml_ploidy), heat = sum(dml_heat), shared = sum(dml_shared)), keyby = lgi]
tab[, all3 := ploidy + heat + shared]
tab[, pct_all3 := round(100 * all3 / sum(all3), 2)][, pct_ploidy_heat := round(100 * (ploidy + heat) / sum(ploidy + heat), 2)]
tab[, pct_all3_placed := round(100 * all3 / sum(all3[!is.na(lgi)]), 2)]
cat("\n== DML per linkage group (lgi = index of NC_047559.1..NC_047568.1; NA = unplaced or mitochondrial) ==\n")
print(tab)
saveRDS(cpg, file.path(scr, "cpg_features.rds"))
