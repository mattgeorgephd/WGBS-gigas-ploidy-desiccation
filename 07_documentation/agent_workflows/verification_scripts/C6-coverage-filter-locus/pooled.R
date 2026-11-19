suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
t5 <- function(f) file.path(repo, "04_output/05_differential_methylation/tables", f)
libs <- readLines("pass/order.txt")
h <- rbindlist(lapply(libs, function(l) fread(sprintf("pass/hist_%s.tsv", l), header = FALSE, col.names = c("cov", "n"))))[cov >= 2, .(n = sum(n)), by = cov][order(cov)]
pooled <- unname(quantile(rep(h$cov, h$n), 99.9 / 100))
cat("single pooled 99.9th percentile over all 22 libraries (rows with >= 2 reads):", pooled, "\n")
z <- readRDS("cand_matrix.rds"); cov <- z$cov; cand <- z$cand
keep <- rowSums(cov >= pooled) == 0
key <- paste(cand$chr, cand$start)[keep]
n <- sapply(c("ploidy", "heat", "shared"), function(s) { d <- fread(t5(sprintf("dml_%s_legacy.csv", s))); sum(paste(d$chr, d$start) %in% key) })
loc <- cand$chr == "NC_047563.1" & cand$start >= 39633323 & cand$start <= 39635520
cat("with one pooled cutoff: tested CpGs", sum(keep), "; removed", sum(!keep), sprintf("(%.2f%%)", 100 * mean(!keep)), "; DML ploidy/heat/shared", n, "; locus CpGs kept", sum(keep & loc), "\n")
