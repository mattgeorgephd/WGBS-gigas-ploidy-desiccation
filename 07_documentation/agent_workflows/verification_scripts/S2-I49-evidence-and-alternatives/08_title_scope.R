# 08_title_scope.R
# Task (c): is "the group differences ... are carried by five libraries ..." supported for the
# draft's own DML lists?  Which libraries carry the draft's ploidy and heat lists?
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
options(width = 220)
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S2-I49-evidence-and-alternatives"
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv")); id <- lib$library_id
heat <- which(lib$treatment == "heat"); ctrl <- which(lib$treatment == "control")
dip <- which(lib$ploidy == "diploid"); tri <- which(lib$ploidy == "triploid")
dh <- intersect(dip, heat); th <- intersect(tri, heat); dc <- intersect(dip, ctrl); tc <- intersect(tri, ctrl)
high <- which(id %in% c("zr3534_3", "zr3534_4", "zr3534_5", "zr3644_14", "zr3644_15"))
m <- readRDS(file.path(scratch, "counts_legacy.rds")); C <- m$C; N <- m$N
pool <- function(idx) 100 * rowSums(C[, idx, drop = FALSE]) / rowSums(N[, idx, drop = FALSE])
dP <- pool(tri) - pool(dip); dH <- pool(heat) - pool(ctrl)
kP <- abs(dP) > 20; kH <- abs(dH) > 20
cat("draft contrasts, CpGs above 20 points: ploidy", sum(kP), "; cohort", sum(kH), "\n")

cat("\n==== the draft's lists without the five high libraries ====\n")
dP5 <- pool(setdiff(tri, high)) - pool(setdiff(dip, high)); dH5 <- pool(setdiff(heat, high)) - pool(setdiff(ctrl, high))
cat("ploidy: above 20 without the five:", sum(abs(dP5) > 20), "; of the", sum(kP), "of the draft,", sum(kP & abs(dP5) > 20 & sign(dP5) == sign(dP)), "remain with the same sign (",
    round(100 * sum(kP & abs(dP5) > 20 & sign(dP5) == sign(dP)) / sum(kP), 1), "%); median |difference| of the draft's CpGs without the five:", round(median(abs(dP5[kP])), 1), "; correlation of the two differences at the draft's CpGs:", round(cor(dP[kP], dP5[kP]), 3), "\n")
cat("cohort: above 20 without the five:", sum(abs(dH5) > 20), "; of the", sum(kH), "of the draft,", sum(kH & abs(dH5) > 20 & sign(dH5) == sign(dH)), "remain with the same sign (",
    round(100 * sum(kH & abs(dH5) > 20 & sign(dH5) == sign(dH)) / sum(kH), 1), "%); correlation at the draft's CpGs:", round(cor(dH[kH], dH5[kH]), 3), "\n")
# reference: leave out five random libraries with the same group membership (3 heat diploids ... cannot: use 2 or 3 per group at random)
set.seed(9)
ref <- t(sapply(1:200, function(b) {
  out <- c(sample(dh, 3), sample(tc, 2))                       # same groups as the five, random members
  d <- pool(setdiff(tri, out)) - pool(setdiff(dip, out)); h <- pool(setdiff(heat, out)) - pool(setdiff(ctrl, out))
  c(ploidy_kept = sum(kP & abs(d) > 20 & sign(d) == sign(dP)) / sum(kP), cohort_kept = sum(kH & abs(h) > 20 & sign(h) == sign(dH)) / sum(kH), is_the_five = setequal(out, high)) }))
ref <- ref[ref[, "is_the_five"] == 0, ]
cat("reference: leaving out 3 random heat diploids and 2 random control triploids (", nrow(ref), "draws other than the five): share of the draft's ploidy CpGs kept: median", round(median(ref[, 1]), 3),
    " range", round(range(ref[, 1]), 3), "; cohort CpGs kept: median", round(median(ref[, 2]), 3), " range", round(range(ref[, 2]), 3), "\n")
# all 10 x 15 = 150 such choices, exhaustively
allc <- expand.grid(a = seq_len(choose(5, 3)), b = seq_len(choose(6, 2)))
ca <- combn(dh, 3, simplify = FALSE); cb <- combn(tc, 2, simplify = FALSE)
ex <- t(sapply(seq_len(nrow(allc)), function(i) { out <- c(ca[[allc$a[i]]], cb[[allc$b[i]]])
  d <- pool(setdiff(tri, out)) - pool(setdiff(dip, out)); h <- pool(setdiff(heat, out)) - pool(setdiff(ctrl, out))
  c(ploidy_kept = sum(kP & abs(d) > 20 & sign(d) == sign(dP)) / sum(kP), cohort_kept = sum(kH & abs(h) > 20 & sign(h) == sign(dH)) / sum(kH), five = setequal(out, high)) }))
cat("exhaustive (150 choices of 3 heat diploids and 2 control triploids): share of the draft's ploidy CpGs kept: the five", round(ex[ex[, 3] == 1, 1], 3), "; others median", round(median(ex[ex[, 3] == 0, 1]), 3),
    " range", round(range(ex[ex[, 3] == 0, 1]), 3), "; rank of the five (1 = fewest kept):", rank(ex[, 1])[ex[, 3] == 1], "of 150\n")
cat("   cohort CpGs kept: the five", round(ex[ex[, 3] == 1, 2], 3), "; others median", round(median(ex[ex[, 3] == 0, 2]), 3), " range", round(range(ex[ex[, 3] == 0, 2]), 3),
    "; rank of the five:", rank(ex[, 2])[ex[, 3] == 1], "of 150\n")

cat("\n==== one library left out: how many of the draft's CpGs remain above 20 points with the same sign ====\n")
loo <- rbindlist(lapply(1:22, function(j) {
  d <- pool(setdiff(tri, j)) - pool(setdiff(dip, j)); h <- pool(setdiff(heat, j)) - pool(setdiff(ctrl, j))
  data.table(library = id[j], group = paste(lib$ploidy[j], lib$treatment[j]), high = j %in% high,
             ploidy_list_without = sum(abs(d) > 20), ploidy_kept = sum(kP & abs(d) > 20 & sign(d) == sign(dP)),
             cohort_list_without = sum(abs(h) > 20), cohort_kept = sum(kH & abs(h) > 20 & sign(h) == sign(dH))) }))
loo[, ploidy_kept_pct := round(100 * ploidy_kept / sum(kP), 1)][, cohort_kept_pct := round(100 * cohort_kept / sum(kH), 1)]
print(loo)

cat("\n==== which contrast does the draft's ploidy list reflect? ====\n")
dHeat <- pool(th) - pool(dh); dCtrl <- pool(tc) - pool(dc)
d25 <- pool(th) - pool(setdiff(dh, high))      # D11, D12 against the five triploids
d35 <- pool(th) - pool(intersect(dh, high))    # the three high diploids against the five triploids
cat("of the", sum(kP), "ploidy CpGs of the draft: above 20 with the same sign in the heat cohort", sum(kP & abs(dHeat) > 20 & sign(dHeat) == sign(dP)),
    "; in the control cohort", sum(kP & abs(dCtrl) > 20 & sign(dCtrl) == sign(dP)),
    "; D11 and D12 against the five heat triploids", sum(kP & abs(d25) > 20 & sign(d25) == sign(dP)),
    "; the three high heat diploids against the five heat triploids", sum(kP & abs(d35) > 20 & sign(d35) == sign(dP)), "\n")
cat("correlation across the draft's ploidy CpGs of the pooled difference with: heat-cohort difference", round(cor(dP[kP], dHeat[kP]), 3), "; control-cohort difference", round(cor(dP[kP], dCtrl[kP]), 3),
    "; D11, D12 against triploids", round(cor(dP[kP], d25[kP]), 3), "; three high diploids against triploids", round(cor(dP[kP], d35[kP]), 3), "\n")
