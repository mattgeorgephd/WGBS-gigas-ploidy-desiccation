suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"; t5 <- file.path(repo, "04_output/05_differential_methylation/tables")
ul <- fread(file.path(t5, "united_cpg_legacy.csv.gz")); hs <- fread("hist_summary.csv")
sf <- max(hs$med10) / hs$med10
cl <- as.matrix(ul[, paste0("coverage", 1:22), with = FALSE])
for (nm in c("cutoff", "cutoff_ge1", "cutoff_ge10")) {
  lim <- round(sf * (ceiling(hs[[nm]]) - 1)); pass <- rowSums(sweep(cl, 2, lim, ">")) == 0
  dp <- fread(file.path(t5, "dml_ploidy_legacy.csv")); dh <- fread(file.path(t5, "dml_heat_legacy.csv"))
  cat(sprintf("percentile base %-12s: CpGs in all 22 libraries at 10x below the cutoff: %d ; legacy ploidy DML at surviving CpGs: %d ; heat: %d ; shared: %d\n", nm, sum(pass),
              sum(pass[dp$united_row]), sum(pass[dh$united_row]), sum(pass[intersect(dp$united_row, dh$united_row)])))
}
