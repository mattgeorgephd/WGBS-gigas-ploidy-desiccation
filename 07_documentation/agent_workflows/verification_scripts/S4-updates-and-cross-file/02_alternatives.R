# S4 check, part 2: tests of the inferences drawn in the updates of I-02, I-04, I-07, I-08.
# Uses the matrices saved by 01_recompute_united.R (legacy united table, N = numCs + numTs).
suppressPackageStartupMessages(library(data.table)); setDTthreads(1)
out <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S4-updates-and-cross-file"
load(file.path(out, "united_legacy_matrices.RData"))
DH <- 1:5; TH <- 6:10; DC <- 11:16; TC <- 17:22; high <- c(3, 4, 5, 18, 19); acute <- c(4, 5, 9, 10)
tol <- 1e-9
rng <- apply(P, 1, function(x) max(x) - min(x))   # a pooled group value lies inside the range of its libraries
keep <- which(rng > 10)                            # a difference of differences above 20 needs a range above 10
cat("CpGs with a library range above 10 points:", length(keep), "; above 20:", sum(rng > 20), "\n")
Ck <- C[keep, ]; Nk <- N[keep, ]
pk <- function(S) { if (length(S) == 1) 100 * Ck[, S] / Nk[, S] else 100 * rowSums(Ck[, S]) / rowSums(Nk[, S]) }
nab <- function(d) sum(abs(d) > 20 + tol)

cat("\n== 1. interaction contrast without the five high libraries, against every within-cohort regrouping of the 17\n")
heat7 <- setdiff(1:10, high); ctrl10 <- setdiff(11:22, high)
ch <- combn(7, 2); cc <- combn(10, 6)          # 2 'diploid' labels among 7 heat, 6 among 10 controls
obs_h <- which(apply(ch, 2, function(s) identical(heat7[s], c(1L, 2L)))); obs_c <- which(apply(cc, 2, function(s) identical(ctrl10[s], 11:16)))
stopifnot(length(obs_h) == 1, length(obs_c) == 1)
dh <- vapply(seq_len(ncol(ch)), function(i) pk(heat7[-ch[, i]]) - pk(heat7[ch[, i]]), numeric(length(keep)))   # 'triploid' minus 'diploid'
dc <- vapply(seq_len(ncol(cc)), function(i) pk(ctrl10[-cc[, i]]) - pk(ctrl10[cc[, i]]), numeric(length(keep)))
res <- matrix(NA_real_, ncol(ch), ncol(cc)); neg <- res; main <- res
for (i in seq_len(ncol(ch))) { for (j in seq_len(ncol(cc))) { d <- dh[, i] - dc[, j]; a <- abs(d) > 20 + tol; res[i, j] <- sum(a); neg[i, j] <- sum(d[a] < 0) } }
cat(sprintf("observed (17 libraries): %d CpGs above 20 points, %.1f%% negative\n", res[obs_h, obs_c], 100 * neg[obs_h, obs_c] / res[obs_h, obs_c]))
cat(sprintf("all %d regroupings: median %.1f ; 5th to 95th %.0f to %.0f ; min %d max %d ; at or above the observed: %d\n",
            length(res), median(res), quantile(res, 0.05), quantile(res, 0.95), min(res), max(res), sum(res >= res[obs_h, obs_c])))
maj <- pmax(neg, res - neg) / res
cat(sprintf("share in the majority direction: observed %.3f ; regroupings median %.3f ; 95th %.3f ; at or above the observed: %d\n",
            maj[obs_h, obs_c], median(maj), quantile(maj, 0.95), sum(maj >= maj[obs_h, obs_c])))
cat("by heat-cohort grouping (the 21 rows), number of control regroupings (of 210) with a count at or above the observed:\n")
print(setNames(rowSums(res >= res[obs_h, obs_c]), apply(ch, 2, function(s) paste(lib$library_id[heat7[s]], collapse = "+"))))
cat("with the heat cohort at its real grouping: control regroupings at or above:", sum(res[obs_h, ] >= res[obs_h, obs_c]), "of 210 ; median", median(res[obs_h, ]), "\n")
cat("with the controls at their real grouping: heat regroupings at or above:", sum(res[, obs_c] >= res[obs_h, obs_c]), "of 21 ; median", median(res[, obs_c]), "\n")
rm(dh, dc); invisible(gc())

cat("\n== 2. acute heat shock inside the heat cohort\n")
c5 <- combn(5, 2)
gr <- expand.grid(i = seq_len(ncol(c5)), j = seq_len(ncol(c5)))
main <- inter <- numeric(nrow(gr)); isobs <- logical(nrow(gr))
for (g in seq_len(nrow(gr))) {
  da <- DH[c5[, gr$i[g]]]; dn <- setdiff(DH, da); ta <- TH[c5[, gr$j[g]]]; tn <- setdiff(TH, ta)
  main[g] <- nab(pk(c(da, ta)) - pk(c(dn, tn)))
  inter[g] <- nab((pk(da) - pk(dn)) - (pk(ta) - pk(tn)))
  isobs[g] <- identical(da, c(4L, 5L)) && identical(ta, c(9L, 10L))
}
cat(sprintf("main contrast (4 with minus 6 without): observed %d ; median of the 100 groupings %.1f ; at or above: %d\n", main[isobs], median(main), sum(main >= main[isobs])))
cat(sprintf("acute by ploidy contrast ((D with - D without) - (T with - T without)): observed %d ; median %.1f ; 5th to 95th %.0f to %.0f ; max %d ; at or above: %d of 100\n",
            inter[isobs], median(inter), quantile(inter, .05), quantile(inter, .95), max(inter), sum(inter >= inter[isobs])))
dd <- vapply(seq_len(ncol(c5)), function(i) nab(pk(DH[c5[, i]]) - pk(DH[-c5[, i]])), numeric(1))
tt <- vapply(seq_len(ncol(c5)), function(i) nab(pk(TH[c5[, i]]) - pk(TH[-c5[, i]])), numeric(1))
names(dd) <- apply(c5, 2, function(s) paste(lib$animal_id[DH[s]], collapse = "+")); names(tt) <- apply(c5, 2, function(s) paste(lib$animal_id[TH[s]], collapse = "+"))
cat("within the diploids, two against three, CpGs above 20 points by pair set apart:\n"); print(sort(dd, decreasing = TRUE))
cat("within the triploids:\n"); print(sort(tt, decreasing = TRUE))
# direction of the acute contrast within each ploidy
dD <- pk(c(4, 5)) - pk(c(1, 2, 3)); dT <- pk(c(9, 10)) - pk(c(6, 7, 8))
cat(sprintf("diploids, D19+D20 minus D11..D13: above 20: %d (%d positive) ; triploids, T19+T20 minus T11..T13: %d (%d positive)\n", nab(dD), sum(dD > 20 + tol), nab(dT), sum(dT > 20 + tol)))
# ploidy inside the heat cohort without the four: 3 against 3
c6 <- combn(6, 3); six <- c(1, 2, 3, 6, 7, 8)
x33 <- vapply(seq_len(ncol(c6)), function(i) nab(pk(six[-c6[, i]]) - pk(six[c6[, i]])), numeric(1))
names(x33) <- apply(c6, 2, function(s) paste(lib$animal_id[six[s]], collapse = "+"))
cat("heat cohort without the four acute libraries, three against three (20 labellings, 10 distinct):\n"); print(sort(x33, decreasing = TRUE))
d33 <- pk(c(6, 7, 8)) - pk(c(1, 2, 3)); cat(sprintf("  observed D11-D13 against T11-T13: %d above 20, %d with diploids higher\n", nab(d33), sum(d33 < -20 - tol)))

cat("\n== 3. how much of a pooled list lies above 20 points inside each cohort, real labels against relabelled\n")
Dall <- function(dip) { tri <- setdiff(1:22, dip); pk(tri) - pk(dip) }
one <- function(hd, cd) {   # hd: 5 heat libraries labelled diploid ; cd: 6 control libraries labelled diploid
  d <- Dall(c(hd, cd)); a <- abs(d) > 20 + tol
  d1 <- pk(setdiff(1:10, hd)) - pk(hd); d2 <- pk(setdiff(11:22, cd)) - pk(cd)
  c(n = sum(a), in_heat = sum(a & abs(d1) > 20 + tol), in_ctrl = sum(a & abs(d2) > 20 + tol), n_heat = sum(abs(d1) > 20 + tol), n_ctrl = sum(abs(d2) > 20 + tol))
}
o <- one(1:5, 11:16); print(o); cat(sprintf("observed: %.1f%% of the pooled list above 20 in the heat cohort, %.1f%% in the control cohort\n", 100 * o["in_heat"] / o["n"], 100 * o["in_ctrl"] / o["n"]))
set.seed(4401)
R <- 300
rr <- t(vapply(seq_len(R), function(r) one(sort(sample(1:10, 5)), sort(sample(11:22, 6))), numeric(5)))
sh <- rr[, "in_heat"] / rr[, "n"]; sc <- rr[, "in_ctrl"] / rr[, "n"]
cat(sprintf("%d random relabellings within cohort: pooled list median %.0f ; share above 20 in the heat cohort: median %.3f (5th to 95th %.3f to %.3f) ; in the control cohort: median %.3f (%.3f to %.3f)\n",
            R, median(rr[, "n"]), median(sh), quantile(sh, .05), quantile(sh, .95), median(sc), quantile(sc, .05), quantile(sc, .95)))
cat(sprintf("relabellings with a heat-cohort share at least the observed %.3f: %d of %d ; with a control-cohort share at most the observed %.3f: %d\n",
            o["in_heat"] / o["n"], sum(sh >= o["in_heat"] / o["n"]), R, o["in_ctrl"] / o["n"], sum(sc <= o["in_ctrl"] / o["n"])))
fwrite(as.data.table(rr), file.path(out, "02_pooled_list_within_cohort_relabelled.csv"))
cat("done\n")
