suppressPackageStartupMessages({library(data.table); library(splines); library(sandwich)})
scr <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C5-feature-odds"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
cpg <- readRDS(file.path(scr, "cpg_features.rds"))
u <- fread(file.path(repo, "04_output/05_differential_methylation/tables/united_cpg_legacy.csv.gz"))
nc <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]); nt <- as.matrix(u[, paste0("numTs", 1:22), with = FALSE])
pm <- 100 * nc / (nc + nt); rm(u, nc, nt); invisible(gc())
rv <- function(m) { k <- ncol(m); (rowSums(m^2) - k * rowMeans(m)^2) / (k - 1) }
pool <- function(g) { v <- sapply(g, function(i) rv(pm[, i, drop = FALSE])); df <- sapply(g, length) - 1; sqrt(as.vector(v %*% df) / sum(df)) }
cpg[, sd_ctrl := pool(list(11:16, 17:22))]      # control cohort only (never heated)
cpg[, sd_heatcoh := pool(list(1:5, 6:10))]      # heat cohort only
e <- cpg[elig_p == TRUE | dml_ploidy | dml_heat]
e[, win := paste(chr, start %/% 100000L)]
for (sdv in c("sd_ctrl", "sd_heatcoh")) {
  cat("\nwithin-group SD from", sdv, "\n")
  e[, sdx := get(sdv)]
  for (s in c("dml_ploidy", "dml_heat")) for (f in c("gene", "CDS", "intron", "TE", "intergenic")) {
    d <- e[elig_p == TRUE | get(s)]; d[, yy := as.integer(get(s))][, ff := as.integer(get(paste0("in_", f)))]
    m <- glm(yy ~ ff + ns(mean_ct, 6) + ns(sdx, 4), binomial, d); b <- coef(m)["ff"]; sc <- sqrt(diag(vcovCL(m, cluster = d$win, type = "HC0"))["ff"])
    cat(sprintf("%-10s %-10s +meth +SD: %.2f (%.2f-%.2f)\n", s, f, exp(b), exp(b - 1.96 * sc), exp(b + 1.96 * sc)))
  }
}
