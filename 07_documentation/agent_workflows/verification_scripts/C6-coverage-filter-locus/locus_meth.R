suppressPackageStartupMessages(library(data.table))
options(width = 230)
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
t5 <- function(f) file.path(repo, "04_output/05_differential_methylation/tables", f)
z <- readRDS("locus.rds"); x <- z$x
pl <- fread("locus_per_library.csv")
L0 <- 39633323L; L1 <- 39635520L
dp <- fread(t5("diffmeth_ploidy_legacy.csv.gz"))[chr == "NC_047563.1" & start >= L0 & start <= L1]
dh <- fread(t5("diffmeth_heat_legacy.csv.gz"))[chr == "NC_047563.1" & start >= L0 & start <= L1]
cat("legacy results at the 54 CpGs: ploidy meth.diff range", round(range(dp$meth.diff), 1), "; DML", sum(dp$qvalue < 0.01 & abs(dp$meth.diff) > 20),
    "(hyper", sum(dp$qvalue < 0.01 & dp$meth.diff > 20), ", hypo", sum(dp$qvalue < 0.01 & dp$meth.diff < -20), "); max q", signif(max(dp$qvalue), 3), "\n")
cat("                               heat   meth.diff range", round(range(dh$meth.diff), 1), "; DML", sum(dh$qvalue < 0.01 & abs(dh$meth.diff) > 20),
    "(hyper", sum(dh$qvalue < 0.01 & dh$meth.diff > 20), ", hypo", sum(dh$qvalue < 0.01 & dh$meth.diff < -20), "); max q", signif(max(dh$qvalue), 3), "\n")
# animal-level methylation at the locus
a <- pl[, .(library_id, ploidy, treatment, hs = acute_heat_shock, state = fifelse(n_over > 0, "high", "low"), ratio = ratio_to_median_all, meth = pooled_meth, median_cov)]
print(a[order(ploidy, treatment, meth)])
cat("\nanimal-level locus methylation (%), unweighted: by ploidy\n"); print(a[, .(n = .N, mean = round(mean(meth), 1), median = median(meth), n_ge30 = sum(meth >= 30)), by = ploidy])
cat("by treatment\n"); print(a[, .(n = .N, mean = round(mean(meth), 1), median = median(meth), n_ge30 = sum(meth >= 30)), by = treatment])
cat("by group\n"); print(a[, .(n = .N, mean = round(mean(meth), 1), median = median(meth), n_ge30 = sum(meth >= 30), n_high_cov = sum(state == "high")), by = .(treatment, ploidy)])
cat("Wilcoxon, animal-level locus methylation: ploidy p =", signif(wilcox.test(meth ~ ploidy, a, exact = FALSE)$p.value, 3),
    "; treatment p =", signif(wilcox.test(meth ~ treatment, a, exact = FALSE)$p.value, 3), "\n")
cat("Fisher, animals with >= 30% methylation by ploidy: p =", signif(fisher.test(table(a$ploidy, a$meth >= 30))$p.value, 3), "\n")
cat("high-coverage state by treatment:\n"); print(table(a$treatment, a$state)); cat("Fisher p =", signif(fisher.test(table(a$treatment, a$state))$p.value, 3), "\n")
cat("high-coverage state by ploidy:\n"); print(table(a$ploidy, a$state)); cat("Fisher p =", signif(fisher.test(table(a$ploidy, a$state))$p.value, 3), "\n")
cat("Wilcoxon, coverage ratio: ploidy p =", signif(wilcox.test(ratio ~ ploidy, a, exact = FALSE)$p.value, 3), "; treatment p =", signif(wilcox.test(ratio ~ treatment, a, exact = FALSE)$p.value, 3), "\n")
# pooled (coverage-weighted) group methylation, as methylKit's wmean, raw counts (no scaling), with and without single libraries
pool <- function(d) d[, 100 * sum(nC) / sum(cov)]
cat("\npooled (coverage-weighted) methylation over the 54 CpGs, raw counts:\n")
cat(" diploid", round(pool(x[ploidy == "diploid"]), 1), " triploid", round(pool(x[ploidy == "triploid"]), 1), " difference (T - D)", round(pool(x[ploidy == "triploid"]) - pool(x[ploidy == "diploid"]), 1), "\n")
cat(" heat", round(pool(x[treatment == "heat"]), 1), " control", round(pool(x[treatment == "control"]), 1), " difference (H - C)", round(pool(x[treatment == "heat"]) - pool(x[treatment == "control"]), 1), "\n")
cat(" share of the heat group's reads at the locus that come from zr3534_4:", round(100 * x[library_id == "zr3534_4", sum(cov)] / x[treatment == "heat", sum(cov)], 1), "%\n")
cat(" without zr3534_4: heat", round(pool(x[treatment == "heat" & library_id != "zr3534_4"]), 1), " difference (H - C)",
    round(pool(x[treatment == "heat" & library_id != "zr3534_4"]) - pool(x[treatment == "control"]), 1),
    "; diploid", round(pool(x[ploidy == "diploid" & library_id != "zr3534_4"]), 1), " difference (T - D)",
    round(pool(x[ploidy == "triploid"]) - pool(x[ploidy == "diploid" & library_id != "zr3534_4"]), 1), "\n")
cat(" unweighted mean of animals: heat", round(a[treatment == "heat", mean(meth)], 1), " control", round(a[treatment == "control", mean(meth)], 1),
    "; diploid", round(a[ploidy == "diploid", mean(meth)], 1), " triploid", round(a[ploidy == "triploid", mean(meth)], 1), "\n")
# per CpG: sign consistency and animal-level tests
x[, pm := 100 * nC / cov]
pc <- x[, .(w_p = wilcox.test(pm ~ ploidy, exact = FALSE)$p.value, w_h = wilcox.test(pm ~ treatment, exact = FALSE)$p.value,
            dmean_p = mean(pm[ploidy == "triploid"]) - mean(pm[ploidy == "diploid"]), dmean_h = mean(pm[treatment == "heat"]) - mean(pm[treatment == "control"]),
            dpool_h_no4 = 100 * sum(nC[treatment == "heat" & library_id != "zr3534_4"]) / sum(cov[treatment == "heat" & library_id != "zr3534_4"]) -
                          100 * sum(nC[treatment == "control"]) / sum(cov[treatment == "control"])), by = start]
cat("\nper CpG (54), animal-level Wilcoxon on percent methylation: ploidy p < 0.05 at", sum(pc$w_p < 0.05), "CpGs (p < 0.01 at", sum(pc$w_p < 0.01), "); treatment p < 0.05 at", sum(pc$w_h < 0.05), "\n")
cat("unweighted mean difference T - D: range", round(range(pc$dmean_p), 1), "; |diff| > 20 at", sum(abs(pc$dmean_p) > 20), "CpGs\n")
cat("unweighted mean difference H - C: range", round(range(pc$dmean_h), 1), "; |diff| > 20 at", sum(abs(pc$dmean_h) > 20), "CpGs\n")
cat("pooled H - C without zr3534_4: range", round(range(pc$dpool_h_no4), 1), "; |diff| > 20 at", sum(abs(pc$dpool_h_no4) > 20), "CpGs\n")
