# S4 check, part 3: "one group apart" (update of I-08; I-49 item 7; I-06 item 6).
# (i) the counts of corrected_dml_group_apart.csv recomputed from corrected_dml_library_values.csv and from the united table;
# (ii) genome-wide: at how many tested CpGs does a set of five heat-cohort libraries lie entirely outside the range of
#      the other 17 libraries, for the real triploid and diploid groups against all 252 sets of five heat libraries?
suppressPackageStartupMessages(library(data.table)); setDTthreads(1)
root <- "/home/user/WGBS-gigas-ploidy-desiccation"
out <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S4-updates-and-cross-file"
load(file.path(out, "united_legacy_matrices.RData"))
v <- fread(file.path(root, "06_diagnostics/library_methylation_level/corrected_dml_library_values.csv"))
cat("rows by mode:\n"); print(v[, .N, by = mode])
vl <- v[mode == "legacy"]
cat("legacy: lists column values (first 6):\n"); print(head(unique(vl$lists)))
inlist <- function(x, l) grepl(l, x, fixed = TRUE)
for (l in c("A1_ploidy_mn_f", "A1_heat_mn_f", "A2_interaction_mn_f", "A4_ploidy_in_control_mn_f", "A4_ploidy_in_heat_mn_f")) {
  s <- vl[inlist(lists, l)]
  cat(sprintf("%-28s n = %3d | ", l, nrow(s))); print(table(s$group_apart_from_all_other_libraries, useNA = "ifany"))
}
# recompute 'apart' for the 125 from the library values in the table itself and from the united table
h <- vl[inlist(lists, "A4_ploidy_in_heat_mn_f")]
libs <- lib$library_id; M <- as.matrix(h[, ..libs])
grp <- list(DH = 1:5, TH = 6:10, DC = 11:16, TC = 17:22)
apart <- function(M, S) { a <- M[, S, drop = FALSE]; b <- M[, -S, drop = FALSE]; (do.call(pmin, as.data.frame(a)) > do.call(pmax, as.data.frame(b))) | (do.call(pmax, as.data.frame(a)) < do.call(pmin, as.data.frame(b))) }
cat("recomputed from the table, the 125 CpGs:", paste(names(grp), vapply(grp, function(S) sum(apart(M, S)), 1L), collapse = " ; "), "\n")
k <- match(paste(h$chr, h$start), paste(chr, pos)); stopifnot(!anyNA(k))
cat("largest difference between table values and percentages recomputed from the united table:", max(abs(M - P[k, ])), "\n")
th_ap <- apart(M, 6:10); dh_ap <- apart(M, 1:5)
gapfun <- function(M, S) { a <- M[, S, drop = FALSE]; b <- M[, -S, drop = FALSE]; pmax(do.call(pmin, as.data.frame(a)) - do.call(pmax, as.data.frame(b)), do.call(pmin, as.data.frame(b)) - do.call(pmax, as.data.frame(a))) }
cat("gap (points) between the group and the nearest other library, triploid-apart CpGs:", paste(round(sort(gapfun(M, 6:10)[th_ap]), 1), collapse = " "), "\n")
cat("gap, diploid-apart CpGs:", paste(round(sort(gapfun(M, 1:5)[dh_ap]), 1), collapse = " "), "\n")
cat("positions of the 17 (chr:start, group):\n"); print(h[th_ap | dh_ap, .(chr, start, group_apart_from_all_other_libraries, range_diploid_heat, range_triploid_heat, range_diploid_control, range_triploid_control)])
# how many of the 125 lie within 1 kb of another of the 125
o <- order(h$chr, h$start); hc <- h$chr[o]; hp <- h$start[o]
cl <- cumsum(c(TRUE, hc[-1] != hc[-length(hc)] | diff(hp) > 1000)); cat("the 125 CpGs form", max(cl), "clusters at 1 kb; the 17 apart form", length(unique(cl[(th_ap | dh_ap)[o]])), "\n")

cat("\n== genome-wide: sets of five heat-cohort libraries entirely outside the range of the other 17 libraries\n")
heat <- 1:10
c5 <- combn(10, 5)
key <- function(s) sum(2^(s - 1))
keys <- apply(c5, 2, key)
ordM <- t(apply(P, 1, order, decreasing = TRUE))     # library indices from highest to lowest
srt <- t(apply(P, 1, sort, decreasing = TRUE))
count_sets <- function(g) {
  top <- srt[, 5] - srt[, 6] > g; bot <- srt[, 17] - srt[, 18] > g
  kt <- rowSums(2^(ordM[top, 1:5, drop = FALSE] - 1)); kb <- rowSums(2^(ordM[bot, 18:22, drop = FALSE] - 1))
  tt <- table(factor(kt, levels = keys)); tb <- table(factor(kb, levels = keys))
  list(top = as.integer(tt), bottom = as.integer(tb), n_top_any = sum(top), n_bot_any = sum(bot), kt = kt, kb = kb)
}
kTH <- key(6:10); kDH <- key(1:5)
for (g in c(0, 10, 20)) {
  r <- count_sets(g); tot <- r$top + r$bottom
  iTH <- which(keys == kTH); iDH <- which(keys == kDH)
  cat(sprintf("gap > %2d points: CpGs with any five libraries apart at the top %d, at the bottom %d ; of these a set of five heat libraries: %d and %d\n", g, r$n_top_any, r$n_bot_any, sum(r$top), sum(r$bottom)))
  cat(sprintf("   triploid heat: %d (higher %d, lower %d) ; diploid heat: %d (higher %d, lower %d)\n", tot[iTH], r$top[iTH], r$bottom[iTH], tot[iDH], r$top[iDH], r$bottom[iDH]))
  oth <- tot[-c(iTH, iDH)]
  cat(sprintf("   the other 250 sets of five heat libraries: median %.1f ; 95th percentile %.1f ; largest %d ; sets with at least the triploid count: %d ; with at least the diploid count: %d\n",
              median(oth), quantile(oth, 0.95), max(oth), sum(oth >= tot[iTH]), sum(oth >= tot[iDH])))
  top5 <- order(tot, decreasing = TRUE)[1:6]
  cat("   six largest sets:", paste(sprintf("%s=%d", apply(c5[, top5], 2, function(s) paste(sub("zr3534_", "", lib$library_id[s]), collapse = ",")), tot[top5]), collapse = " | "), "\n")
  if (g == 0) {
    # all sets of five of the 22 libraries, not only heat-cohort sets
    allk <- c(r$kt, r$kb); ta <- sort(table(allk), decreasing = TRUE)
    cat(sprintf("   over all sets of five of the 22 libraries: %d distinct sets occur among 26,334 ; counts of the ten most frequent: %s ; rank of triploid heat: %d ; of diploid heat: %d\n",
                length(ta), paste(as.integer(ta[1:10]), collapse = " "), which(names(ta) == as.character(kTH)), which(names(ta) == as.character(kDH))))
    dec <- function(kk) paste(lib$library_id[which(bitwAnd(as.integer(kk), 2^(0:21)) > 0)], collapse = ",")
    cat("   the five most frequent sets:\n"); for (z in 1:5) cat("     ", as.integer(ta[z]), ":", dec(names(ta)[z]), "\n")
  }
}
# where are the triploid-apart CpGs, and how clustered?
r <- count_sets(0)
top <- srt[, 5] - srt[, 6] > 0; bot <- srt[, 17] - srt[, 18] > 0
isTH <- rep(FALSE, nrow(P)); isTH[which(top)[r$kt == kTH]] <- TRUE; isTH[which(bot)[r$kb == kTH]] <- TRUE
isDH <- rep(FALSE, nrow(P)); isDH[which(top)[r$kt == kDH]] <- TRUE; isDH[which(bot)[r$kb == kDH]] <- TRUE
clus <- function(sel) { o <- order(chr[sel], pos[sel]); a <- chr[sel][o]; b <- pos[sel][o]; sum(c(TRUE, a[-1] != a[-length(a)] | diff(b) > 1000)) }
cat(sprintf("triploid heat apart at %d CpGs in %d clusters (1 kb) ; diploid heat apart at %d CpGs in %d clusters\n", sum(isTH), clus(isTH), sum(isDH), clus(isDH)))
cat("by chromosome, triploid heat apart:\n"); print(table(chr[isTH]))
cat("by chromosome, diploid heat apart:\n"); print(table(chr[isDH]))
gTH <- pmax(apply(P[isTH, 6:10, drop = FALSE], 1, min) - apply(P[isTH, -(6:10), drop = FALSE], 1, max), apply(P[isTH, -(6:10), drop = FALSE], 1, min) - apply(P[isTH, 6:10, drop = FALSE], 1, max))
cat("gap quartiles, triploid-apart CpGs:", paste(round(quantile(gTH, c(0, .25, .5, .75, 1)), 1), collapse = " "), "\n")
cat("done\n")
