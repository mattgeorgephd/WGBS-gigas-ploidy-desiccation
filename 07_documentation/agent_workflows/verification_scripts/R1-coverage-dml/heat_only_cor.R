# Side question (I-16): can the manuscript's "0.87 to 0.93" be a correlation range of a subset analysis (heat libraries only, as in the
# older legacy script 08_archive/code/2_WGBS_Methylkit.R)? One coverage file in memory at a time.
suppressPackageStartupMessages(library(data.table)); setDTthreads(2)
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
source(file.path(repo, "02_R_functions/coverage_io.R"))
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv")); setorder(lib, legacy_sample_index)
set_name <- commandArgs(TRUE)[1]; idx <- if (set_name == "heat") which(lib$treatment == "heat") else which(lib$treatment == "control")
lev <- sort(fread("union/union_by_seq.tsv", header = FALSE)$V1)
files <- file.path("/home/user/gannet_cache/cov_files", lib$cov_file[idx]); n <- length(idx)
common <- NULL
for (j in seq_len(n)) { d <- read_cov(files[j])[-1L]; k <- cpg_key(d$chr, d$start, lev); common <- if (is.null(common)) k else common[common %in% k]; rm(d, k); gc(verbose = FALSE)
  cat(sprintf("pass 1 %s: common %d\n", lib$library_id[idx[j]], length(common))) }
cv <- matrix(0L, length(common), n); nc <- matrix(0L, length(common), n); med <- list()
for (j in seq_len(n)) { d <- read_cov(files[j])[-1L]; d[, coverage := nC + nT]; m <- match(common, cpg_key(d$chr, d$start, lev)); cv[, j] <- d$coverage[m]; nc[, j] <- d$nC[m]
  med[[j]] <- vapply(c(1, 3, 5, 10, 20), function(x) median(d$coverage[d$coverage >= x]), numeric(1)); rm(d, m); gc(verbose = FALSE) }
med <- do.call(cbind, med); rownames(med) <- c(1, 3, 5, 10, 20)
cat(sprintf("%s-only set: %d libraries, %d CpGs with >= 1 read in all of them\n", set_name, n, length(common)))
for (x in c(1, 3, 5, 10, 20)) {
  rows <- which(rowSums(cv >= x) == n); sf <- max(med[as.character(x), ]) / med[as.character(x), ]
  pm <- vapply(seq_len(n), function(j) { a <- cv[rows, j]; b <- nc[rows, j]; nw <- round(sf[j] * a); cs <- round(nw * (b / a)); ts <- round(nw * ((a - b) / a)); cs / (cs + ts) }, numeric(length(rows)))
  r <- cor(pm, use = "na.or.complete"); off <- r[upper.tri(r)]
  raw <- cor(nc[rows, ] / cv[rows, ]); offr <- raw[upper.tri(raw)]
  cat(sprintf("  %2dx: %8d CpGs ; Pearson r (methylKit-normalized) %.4f to %.4f ; on raw counts %.4f to %.4f\n", x, length(rows), min(off), max(off), min(offr), max(offr)))
}
