suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
tdir <- file.path(repo, "04_output/05_differential_methylation/tables")
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv")); setorder(lib, legacy_sample_index)
fl <- fread(file.path(tdir, "library_filtering_legacy.csv"))
idx <- fread("rawcheck/idx.csv")$united_row
u <- fread(file.path(tdir, "united_cpg_legacy.csv.gz"))[idx]
bad <- 0; tot <- 0; badlib <- integer(0)
for (j in 1:22) {
  r <- fread(sprintf("rawcheck/%s.sel", lib$cov_file[j]), header = FALSE, col.names = c("chr", "start", "end", "pct", "nC", "nT"))
  m <- match(paste(u$chr, u$start), paste(r$chr, r$start))
  if (anyNA(m)) { cat("library", j, ": united CpGs absent from the raw file:", sum(is.na(m)), "\n"); next }
  r <- r[m]; cv <- r$nC + r$nT; sf <- max(fl$median_coverage) / fl$median_coverage[j]
  nc <- round(sf * cv); e <- data.table(coverage = nc, numCs = round(nc * (r$nC / cv)), numTs = round(nc * (r$nT / cv)))
  o <- u[, c(paste0("coverage", j), paste0("numCs", j), paste0("numTs", j)), with = FALSE]
  d <- sum(e$coverage != o[[1]]) + sum(e$numCs != o[[2]]) + sum(e$numTs != o[[3]])
  if (d > 0) badlib <- c(badlib, j)
  bad <- bad + d; tot <- tot + 3 * nrow(u)
  if (any(cv < 10)) cat("library", j, ": raw coverage < 10 at a united CpG\n")
  if (any(r$end != u$end)) cat("library", j, ": end differs\n")
}
cat(sprintf("%d united CpGs x 22 libraries: %d of %d scaled counts differ from recomputation out of the raw files (libraries with differences: %s)\n", nrow(u), bad, tot, paste(badlib, collapse = ",")))
# sanity: would a shifted library order be detected? compare column j with raw library j+1
r1 <- fread(sprintf("rawcheck/%s.sel", lib$cov_file[2]), header = FALSE, col.names = c("chr", "start", "end", "pct", "nC", "nT"))
m <- match(paste(u$chr, u$start), paste(r1$chr, r1$start)); cv <- (r1$nC + r1$nT)[m]
cat("control: column 1 against raw library 2 mismatches:", sum(round(max(fl$median_coverage) / fl$median_coverage[2] * cv) != u$coverage1), "of", nrow(u), "\n")
