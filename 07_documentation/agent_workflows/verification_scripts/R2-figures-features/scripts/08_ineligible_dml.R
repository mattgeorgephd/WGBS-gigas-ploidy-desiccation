suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"; t05 <- file.path(repo, "04_output/05_differential_methylation/tables")
u <- fread(file.path(t05, "united_cpg_legacy.csv.gz"))
lib <- read.csv(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv")); lib <- lib[order(lib$legacy_sample_index), ]
for (k in list(c("NC_047560.1", 7406728, "ploidy"), c("NC_047562.1", 38931611, "heat"))) {
  r <- which(u$chr == k[1] & u$start == as.integer(k[2])); d <- fread(file.path(t05, sprintf("dml_%s_legacy.csv", k[3])))[united_row == r]
  nC <- unlist(u[r, paste0("numCs", 1:22), with = FALSE]); nT <- unlist(u[r, paste0("numTs", 1:22), with = FALSE])
  pm <- 100 * nC / (nC + nT); g <- if (k[3] == "ploidy") lib$ploidy else lib$treatment
  cat(sprintf("\n%s:%s (%s DML) row %d  meth.diff %.2f  q %.3g\n", k[1], k[2], k[3], r, d$meth.diff, d$qvalue))
  print(data.frame(lib = lib$library_id, group = g, numCs = nC, numTs = nT, pct = round(pm, 1)), row.names = FALSE)
  cat(sprintf("unweighted mean of library percentages: %.2f; by group: %s\n", mean(pm), paste(names(tapply(pm, g, mean)), round(tapply(pm, g, mean), 2), collapse = ", ")))
  cat(sprintf("pooled (coverage-weighted) percentage by group: %s; difference %.2f\n", paste(names(tapply(nC, g, sum)), round(100 * tapply(nC, g, sum) / (tapply(nC, g, sum) + tapply(nT, g, sum)), 2), collapse = ", "),
              diff(100 * tapply(nC, g, sum) / (tapply(nC, g, sum) + tapply(nT, g, sum)))))
}
# how many analysed CpGs outside the eligible range have |pooled group difference| > 20 for either contrast?
nC <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]); nT <- as.matrix(u[, paste0("numTs", 1:22), with = FALSE])
pm <- 100 * nC / (nC + nT); mp <- rowMeans(pm); elig <- mp >= 10 & mp <= 90
pooled <- function(idx) 100 * rowSums(nC[, idx]) / (rowSums(nC[, idx]) + rowSums(nT[, idx]))
dpl <- pooled(lib$ploidy == "triploid") - pooled(lib$ploidy == "diploid"); dht <- pooled(lib$treatment == "heat") - pooled(lib$treatment == "control")
cat(sprintf("\nlegacy: CpGs with |pooled ploidy difference| > 20: %d, of which outside the eligible range: %d; heat: %d, outside: %d\n", sum(abs(dpl) > 20), sum(abs(dpl) > 20 & !elig), sum(abs(dht) > 20), sum(abs(dht) > 20 & !elig)))
dm <- fread(file.path(t05, "diffmeth_ploidy_legacy.csv.gz")); cat("max |my pooled ploidy diff - methylKit meth.diff|:", max(abs(dpl - dm$meth.diff)), "\n")
dm <- fread(file.path(t05, "diffmeth_heat_legacy.csv.gz")); cat("max |my pooled heat diff - methylKit meth.diff|:", max(abs(dht - dm$meth.diff)), "\n")
