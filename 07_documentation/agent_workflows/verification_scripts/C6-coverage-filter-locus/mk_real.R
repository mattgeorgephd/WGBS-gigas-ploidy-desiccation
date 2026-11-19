# methylKit's own methRead + filterByCoverage on ONE real coverage file, against the pipeline's records and my awk extraction.
suppressPackageStartupMessages({library(methylKit); library(data.table)})
lib <- commandArgs(trailingOnly = TRUE)[1]
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
f <- sprintf("/home/user/gannet_cache/cov_files/%s_R1.CpG_report.merged_CpG_evidence.cov", lib)
t5 <- function(x) file.path(repo, "04_output/05_differential_methylation/tables", x)
recL <- fread(t5("library_filtering_legacy.csv"))[library_id == lib]; recC <- fread(t5("library_filtering_corrected.csv"))[library_id == lib]
mine <- fread("cutoffs_mine.csv")[library_id == lib]
line1 <- fread("first_lines.tsv", header = FALSE)[V1 == lib]
out <- list(library = lib)
mem <- function() round(sum(gc(verbose = FALSE)[, 6]), 0)       # max MB used so far
# ---- A: the legacy call: default header, mincov = 2, then the misnamed upper filter ----
A <- suppressMessages(methRead(f, sample.id = lib, assembly = "Roslin", pipeline = "bismarkCoverage", mincov = 2))
out$A_rows <- nrow(A); out$A_first <- sprintf("%s:%d", A$chr[1], A$start[1])
out$A_has_line1 <- any(A$chr == line1$V2 & A$start == line1$V3)
Af <- filterByCoverage(A, lo.count = 10, lo.perc = NULL, high.count = NULL, high.perc = 99.9)
out$A_legacy_filter_rows <- nrow(Af); out$A_legacy_filter_maxcov <- max(Af$coverage); out$A_legacy_median <- median(Af$coverage)
rm(A, Af); gc(verbose = FALSE)
# ---- B: complete read (header = FALSE), then the filter with the right argument name ----
B <- suppressMessages(methRead(f, sample.id = lib, assembly = "Roslin", pipeline = "bismarkCoverage", mincov = 2, header = FALSE))
out$B_rows <- nrow(B); out$B_first <- sprintf("%s:%d", B$chr[1], B$start[1])
out$B_cutoff <- unname(quantile(B$coverage, 99.9 / 100))
Bf <- filterByCoverage(B, lo.count = 10, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)
out$B_filter_rows <- nrow(Bf); out$B_filter_maxcov <- max(Bf$coverage); out$B_filter_mincov <- min(Bf$coverage); out$B_median <- median(Bf$coverage)
d <- as.data.table(getData(B)[B$coverage >= 10, c("chr", "start", "coverage")]); rm(B); gc(verbose = FALSE)
kept <- as.data.table(getData(Bf)[, c("chr", "start")]); kept[, kept := TRUE]
removed <- kept[d, on = .(chr, start)][is.na(kept), .(chr, start, coverage)]; rm(kept)
out$B_removed <- nrow(removed); out$B_removed_mincov <- min(removed$coverage)
hi <- fread(sprintf("pass/hi_%s.tsv", lib), header = FALSE, col.names = c("chr", "start", "nC", "nT"))[nC + nT >= mine$cutoff_mine]
out$removed_set_equals_awk <- setequal(paste(removed$chr, removed$start), paste(hi$chr, hi$start)) && nrow(removed) == nrow(hi)
rm(d); out$peak_MB <- mem()
cat(sprintf("\n===== %s =====\n", lib))
cat(sprintf("file lines %d (first line %s:%d, coverage %d)\n", mine$n_lines, line1$V2, line1$V3, line1$V6 + line1$V7))
cat(sprintf("A  default header, mincov 2: rows %d (pipeline legacy n_cpg_imported %d; rows with >= 2 reads in file %d) | first row %s | line 1 present: %s\n",
            out$A_rows, recL$n_cpg_imported, mine$n_ge2, out$A_first, out$A_has_line1))
cat(sprintf("A  legacy filter call (high.perc): rows %d (pipeline legacy n_cpg_kept %d) | max coverage kept %d (pipeline %d) | median %g (pipeline %g)\n",
            out$A_legacy_filter_rows, recL$n_cpg_kept, out$A_legacy_filter_maxcov, recL$max_coverage_kept, out$A_legacy_median, recL$median_coverage))
cat(sprintf("B  header = FALSE, mincov 2: rows %d (pipeline corrected n_cpg_imported %d) | first row %s\n", out$B_rows, recC$n_cpg_imported, out$B_first))
cat(sprintf("B  methylKit cutoff %.10f | pipeline %.10f | mine %.10f | identical to pipeline csv value within 1e-9: %s\n", out$B_cutoff, recC$upper_cutoff, mine$cutoff_mine,
            abs(out$B_cutoff - recC$upper_cutoff) < 1e-9))
cat(sprintf("B  filterByCoverage(lo.count = 10, hi.perc = 99.9): rows %d (pipeline n_cpg_kept %d) | coverage kept %d to %d (pipeline max %d) | median %g (pipeline %g)\n",
            out$B_filter_rows, recC$n_cpg_kept, out$B_filter_mincov, out$B_filter_maxcov, recC$max_coverage_kept, out$B_median, recC$median_coverage))
cat(sprintf("B  removed by the upper filter: %d (pipeline n_removed_by_upper %d) | smallest removed coverage %d | same positions as awk rule (cov >= cutoff): %s\n",
            out$B_removed, recC$n_removed_by_upper, out$B_removed_mincov, out$removed_set_equals_awk))
cat(sprintf("ALL AGREE: %s | peak memory %d MB\n",
            out$A_rows == recL$n_cpg_imported && !out$A_has_line1 && out$A_legacy_filter_rows == recL$n_cpg_kept && out$A_legacy_filter_maxcov == recL$max_coverage_kept &&
            out$B_rows == recC$n_cpg_imported && abs(out$B_cutoff - recC$upper_cutoff) < 1e-9 && out$B_filter_rows == recC$n_cpg_kept &&
            out$B_filter_maxcov == recC$max_coverage_kept && out$B_removed == recC$n_removed_by_upper && out$removed_set_equals_awk &&
            out$A_legacy_median == recL$median_coverage && out$B_median == recC$median_coverage, out$peak_MB))
