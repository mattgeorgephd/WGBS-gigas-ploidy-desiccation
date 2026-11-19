suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
t5 <- function(f) file.path(repo, "04_output/05_differential_methylation/tables", f)
ord <- readLines("pass/order.txt")
cut <- fread("cutoffs_mine.csv")
cand <- fread("pass/cand_after_zr3534_8.tsv", header = FALSE, sep = "\t", col.names = c("chr", "start", "s"))
cat("CpGs with >= 10 reads in all 22 libraries (complete read of every file, no upper filter):", nrow(cand), "\n")
parts <- tstrsplit(cand$s, ",", fixed = TRUE); stopifnot(length(parts) == 22)
cov <- matrix(0L, nrow(cand), 22, dimnames = list(NULL, ord)); nC <- cov
for (j in 1:22) { x <- tstrsplit(parts[[j]], ":", fixed = TRUE); cov[, j] <- as.integer(x[[1]]); nC[, j] <- as.integer(x[[2]]) }
rm(parts); cand[, s := NULL]
stopifnot(all(cov >= 10))
first <- fread("first_lines.tsv", header = FALSE, col.names = c("library_id", "chr", "start", "end", "pct", "nC", "nT"))
# legacy: a library's first-line CpG does not count for that library
fl <- merge(cand[, .(chr, start, row = .I)], unique(first[, .(chr, start)]), by = c("chr", "start"))
cat("first-line CpGs (of any file) among the 480,579:", nrow(fl), "\n")
cm <- function(col) cut[match(ord, library_id)][[col]]
keep <- function(cutoffs) rowSums(cov >= matrix(cutoffs, nrow(cov), 22, byrow = TRUE)) == 0
k2 <- keep(cm("cutoff_mine"))
cat("after removing CpGs at or above the 99.9th percentile (base: >= 2 reads) in any library:", sum(k2), " removed:", sum(!k2),
    sprintf("(%.2f%%)", 100 * mean(!k2)), "\n")
cat("  same with strict '>' (remove only above the cutoff):", sum(rowSums(cov > matrix(cm("cutoff_mine"), nrow(cov), 22, byrow = TRUE)) == 0), "\n")
cat("  percentile base >= 1 read (methRead(mincov = 1), as in 1_coverage_analysis.R):", sum(keep(cm("cutoff_base1"))), "\n")
cat("  percentile base >= 10 reads (methRead default mincov = 10):", sum(keep(cm("cutoff_base10"))), "\n")
cat("  cutoff from the legacy import (line 1 lost):", sum(keep(cm("cutoff_line1lost"))), "\n")
# how many libraries exceed, for removed CpGs
nlib_over <- rowSums(cov >= matrix(cm("cutoff_mine"), nrow(cov), 22, byrow = TRUE))
cat("removed CpGs by number of libraries at or above their cutoff:\n"); print(table(nlib_over[nlib_over > 0]))
cat(sprintf("removed CpGs over the cutoff in exactly 1 library: %d (%.1f%% of removed); in all 22: %d (%.1f%%)\n", sum(nlib_over == 1),
            100 * sum(nlib_over == 1) / sum(nlib_over > 0), sum(nlib_over == 22), 100 * sum(nlib_over == 22) / sum(nlib_over > 0)))
# ---- compare with the pipeline's united tables: exact set equality ----
for (mode in c("legacy", "corrected")) {
  u <- fread(t5(sprintf("united_cpg_%s.csv.gz", mode)), select = c("chr", "start"))
  mine <- if (mode == "legacy") cand else cand[k2]
  cat(sprintf("[%s] pipeline united rows %d; mine %d; identical position sets: %s\n", mode, nrow(u), nrow(mine),
              setequal(paste(u$chr, u$start), paste(mine$chr, mine$start))))
}
# ---- raw counts against the united table values (methylKit median scaling) for the corrected mode ----
lf <- fread(t5("library_filtering_corrected.csv"))
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv"))
u <- fread(t5("united_cpg_corrected.csv.gz"))
m <- cand[k2][, row := which(k2)]
u[, urow := .I]; mm <- merge(u[, .(chr, start, urow)], m, by = c("chr", "start"))
ok_all <- TRUE
for (i in seq_len(nrow(lib))) {
  l <- lib$library_id[i]; sf <- lf[library_id == l, scaling_factor]
  cv <- cov[mm$row, l]; c1 <- nC[mm$row, l]
  new_cov <- round(sf * cv); new_c <- round(new_cov * (c1 / cv)); new_t <- round(new_cov * ((cv - c1) / cv))
  ok <- identical(as.numeric(u[[paste0("coverage", i)]][mm$urow]), as.numeric(new_cov)) &&
        identical(as.numeric(u[[paste0("numCs", i)]][mm$urow]), as.numeric(new_c)) &&
        identical(as.numeric(u[[paste0("numTs", i)]][mm$urow]), as.numeric(new_t))
  ok_all <- ok_all && ok
}
cat("corrected united table: coverage, numCs, numTs of all 22 libraries reproduced from the raw files with the table's scaling factors:", ok_all, "\n")
# medians of each library's filtered set, from my histograms -> scaling factors
med <- sapply(lib$library_id, function(l) { h <- fread(sprintf("pass/hist_%s.tsv", l), header = FALSE, col.names = c("cov", "n"))
  ct <- cut[library_id == l, cutoff_mine]; a <- h[cov >= 10]; b <- h[cov >= 10 & cov < ct]
  c(legacy = median(rep(a$cov, a$n)), corrected = median(rep(b$cov, b$n)), kept_corrected = sum(b$n), removed = sum(a$n) - sum(b$n)) })
cat("median coverage of each library's filtered set equals the pipeline record (corrected):", all(med["corrected", ] == lf$median_coverage),
    "; n kept equal:", all(med["kept_corrected", ] == lf$n_cpg_kept), "; n removed equal:", all(med["removed", ] == lf$n_removed_by_upper), "\n")
cat("scaling factors reproduced:", isTRUE(all.equal(max(med["corrected", ]) / med["corrected", ], lf$scaling_factor, check.attributes = FALSE)),
    "; same in both modes:", all(med["legacy", ] == med["corrected", ]), "\n")
saveRDS(list(cand = cand, cov = cov, nC = nC, k2 = k2, nlib_over = nlib_over), "cand_matrix.rds")
