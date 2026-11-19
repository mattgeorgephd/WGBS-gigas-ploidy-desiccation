# Leave-one-library-out at the 54 locus CpGs, with the normalized counts the legacy model used (united_cpg_legacy).
suppressPackageStartupMessages(library(data.table))
options(width = 200)
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
t5 <- function(f) file.path(repo, "04_output/05_differential_methylation/tables", f)
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv"))
trip <- as.integer(lib$ploidy == "triploid"); heat <- as.integer(lib$treatment == "heat")
u <- fread(t5("united_cpg_legacy.csv.gz"))[chr == "NC_047563.1" & start >= 39633323 & start <= 39635520]
C <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]); T <- as.matrix(u[, paste0("numTs", 1:22), with = FALSE])
one <- function(keep, g, cv) t(vapply(seq_len(nrow(C)), function(i) {
  cs <- C[i, keep]; ts <- T[i, keep]; gg <- g[keep]; cc <- cv[keep]
  f1 <- suppressWarnings(glm(cbind(cs, ts) ~ gg + cc, family = binomial)); f0 <- suppressWarnings(glm(cbind(cs, ts) ~ cc, family = binomial))
  c(p = pchisq(f0$deviance - f1$deviance, 1, lower.tail = FALSE),
    diff = 100 * (sum(cs[gg == 1]) / sum((cs + ts)[gg == 1]) - sum(cs[gg == 0]) / sum((cs + ts)[gg == 0]))) }, numeric(2)))
pth <- 2e-4      # p at which q reaches 0.01 in the legacy/corrected tables (ploidy 2.0e-4, heat 4.1e-4); conservative for both
full_p <- one(rep(TRUE, 22), trip, heat); full_h <- one(rep(TRUE, 22), heat, trip)
cat("all 22 libraries: ploidy |diff| > 20 & p <", pth, ":", sum(abs(full_p[, "diff"]) > 20 & full_p[, "p"] < pth), "; heat:", sum(abs(full_h[, "diff"]) > 20 & full_h[, "p"] < pth), "\n")
res <- rbindlist(lapply(1:22, function(j) { k <- rep(TRUE, 22); k[j] <- FALSE; rp <- one(k, trip, heat); rh <- one(k, heat, trip)
  data.table(left_out = lib$library_id[j], ploidy = lib$ploidy[j], treatment = lib$treatment[j],
             share_of_reads_pct = round(100 * sum(C[, j] + T[, j]) / sum(C + T), 1),
             ploidy_dml = sum(abs(rp[, "diff"]) > 20 & rp[, "p"] < pth), ploidy_diff_median = round(median(rp[, "diff"]), 1),
             heat_dml = sum(abs(rh[, "diff"]) > 20 & rh[, "p"] < pth), heat_diff_median = round(median(rh[, "diff"]), 1)) }))
print(res)
cat("share of each group's normalized reads at the locus from its largest library: heat", round(100 * max(colSums(C + T)[heat == 1]) / sum((C + T)[, heat == 1]), 1),
    "% (", lib$library_id[heat == 1][which.max(colSums(C + T)[heat == 1])], "); control", round(100 * max(colSums(C + T)[heat == 0]) / sum((C + T)[, heat == 0]), 1),
    "%; diploid", round(100 * max(colSums(C + T)[trip == 0]) / sum((C + T)[, trip == 0]), 1), "%; triploid", round(100 * max(colSums(C + T)[trip == 1]) / sum((C + T)[, trip == 1]), 1), "%\n")
# only the 9 libraries with ordinary coverage at the locus
pl <- fread("locus_per_library.csv"); low <- lib$library_id %in% pl[n_over == 0, library_id]
rp <- one(low, trip, heat); rh <- one(low, heat, trip)
cat("only the 9 ordinary-coverage libraries (", sum(low & trip == 0), "diploid,", sum(low & trip == 1), "triploid;", sum(low & heat == 1), "heat,", sum(low & heat == 0), "control ):",
    "ploidy |diff| > 20 & p <", pth, ":", sum(abs(rp[, "diff"]) > 20 & rp[, "p"] < pth), " median diff", round(median(rp[, "diff"]), 1),
    "; heat:", sum(abs(rh[, "diff"]) > 20 & rh[, "p"] < pth), " median diff", round(median(rh[, "diff"]), 1), "\n")
