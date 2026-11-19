suppressPackageStartupMessages({library(data.table); library(methylKit)})
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
source(file.path(repo, "02_R_functions/coverage_io.R")); source(file.path(repo, "02_R_functions/methylkit_helpers.R"))
ok <- function(label, cond) cat(sprintf("%-100s %s\n", label, if (isTRUE(cond)) "OK" else "*** FAIL ***"))
# 1. key round trip
lev <- sprintf("seq%03d", 1:400); set.seed(5)
chr <- sample(lev, 1e6, TRUE); st <- c(0L, 999999999L, sample.int(999999999L, 1e6 - 2))
k <- cpg_key(chr, st, lev); p <- key_to_position(k, lev)
ok("cpg_key / key_to_position round trip (400 sequences, starts 0 to 999,999,999)", identical(p$chr, chr) && identical(p$start, st) && !anyDuplicated(k[!duplicated(paste(chr, st))]))
ok("cpg_key stops on an unknown sequence name", inherits(try(cpg_key("zz", 1L, lev), silent = TRUE), "try-error"))
# 2. cutoff vs methylKit::filterByCoverage on synthetic counts, several percentiles
mk <- function(n) { cv <- pmax(1L, as.integer(rnbinom(n, mu = 9, size = 1.2))); nC <- rbinom(n, cv, rbeta(n, 0.3, 1.5)); data.frame(chr = "c", start = seq_len(n) * 3L, end = seq_len(n) * 3L + 2L, strand = "*", coverage = cv, numCs = nC, numTs = cv - nC, stringsAsFactors = FALSE) }
allok <- TRUE
for (pc in c(99.9, 99.5, 99, 95, 90.1)) for (rep in 1:5) {
  d <- mk(20011); imp <- d[d$coverage >= 2, ]
  obj <- new("methylRaw", imp, sample.id = "s", assembly = "a", context = "CpG", resolution = "base")
  nat <- getData(filterByCoverage(obj, lo.count = 10, hi.perc = pc))
  cut <- upper_coverage_cutoff(d$coverage, 2, pc); mine <- d[d$coverage >= 10 & d$coverage < cut, ]
  allok <- allok && identical(nat$start, mine$start)
}
ok("upper_coverage_cutoff + '< cutoff' keeps the rows filterByCoverage(hi.perc) keeps (25 synthetic cases)", allok)
# 3. scaling vs normalizeCoverage, bitwise, including ties
allok <- TRUE; nties <- 0
for (rep in 1:20) {
  objs <- lapply(1:4, function(i) { d <- mk(30000); d$coverage <- d$coverage * sample(1:3, 1); d$numCs <- pmin(d$numCs, d$coverage); d$numTs <- d$coverage - d$numCs; d <- d[d$coverage >= 10, ]
    new("methylRaw", d, sample.id = as.character(i), assembly = "a", context = "CpG", resolution = "base") })
  nrm <- normalizeCoverage(new("methylRawList", objs, treatment = c(0, 1, 0, 1)))
  med <- sapply(objs, function(o) median(o$coverage)); sf <- max(med) / med
  for (i in 1:4) { o <- getData(objs[[i]]); s <- scale_counts_like_methylkit(o$coverage, o$numCs, o$numTs, sf[i]); n <- getData(nrm[[i]])
    allok <- allok && identical(s$coverage, n$coverage) && identical(s$numCs, n$numCs) && identical(s$numTs, n$numTs)
    nties <- nties + sum(abs((sf[i] * o$coverage) %% 1 - 0.5) < 1e-9) }
}
ok(sprintf("scale_counts_like_methylkit identical to normalizeCoverage (80 synthetic libraries, %d coverage ties)", nties), allok)
# 4. read_filtered_counts: first line with >= 10 reads is lost in legacy mode only; cutoff the same in both modes
tf <- tempfile(fileext = ".cov"); d <- mk(5000); d$coverage[1] <- 40L; d$numCs[1] <- 10L; d$numTs[1] <- 30L
fwrite(data.table(d$chr, d$start, d$end, 100 * d$numCs / d$coverage, d$numCs, d$numTs), tf, sep = "\t", col.names = FALSE)
cfg <- list(import_min_coverage = 2, min_coverage = 10, upper_coverage_percentile = 99.9)
a <- read_filtered_counts(tf, "legacy", cfg, character(0)); b <- read_filtered_counts(tf, "corrected", cfg, character(0))
ok("legacy drops the first line, corrected keeps it; same cutoff", !(cpg_key("c", d$start[1], "c") %in% a$counts$pos_key) && (cpg_key("c", d$start[1], "c") %in% b$counts$pos_key) && a$record$upper_cutoff == b$record$upper_cutoff)
ok("unknown mode stops", inherits(try(read_filtered_counts(tf, "Legacy", cfg, character(0)), silent = TRUE), "try-error"))
