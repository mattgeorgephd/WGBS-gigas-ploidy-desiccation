suppressPackageStartupMessages(library(methylKit))
set.seed(1)
S <- "toy"
# ---- toy headerless Bismark coverage files (zero-based merged CpG format) ----
mk <- function(path, n = 5000, seed = 1) {
  set.seed(seed)
  cov <- pmax(1L, rnbinom(n, size = 2, mu = 12))
  cov[sample(n, 8)] <- c(400L, 350L, 300L, 250L, 200L, 180L, 170L, 160L)   # a few extreme sites
  nC <- rbinom(n, cov, 0.3)
  d <- data.frame(chr = "chrA", start = seq(100L, by = 7L, length.out = n), stringsAsFactors = FALSE)
  d$end <- d$start + 2L; d$pct <- sprintf("%.6f", 100 * nC / cov); d$nC <- nC; d$nT <- cov - nC
  write.table(d, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  d
}
f1 <- file.path(S, "a.cov"); f2 <- file.path(S, "b.cov")
d1 <- mk(f1, seed = 1); d2 <- mk(f2, seed = 2)
cat("lines in toy file a:", nrow(d1), " first line:", paste(d1[1, ], collapse = " "), "\n")

# ---- (I-38) methRead default header ----
cat("\nformals(methRead)$header =", deparse(formals(methRead)$header), "\n")
o_def <- suppressMessages(methRead(f1, sample.id = "a", assembly = "x", pipeline = "bismarkCoverage", mincov = 1))
o_nohdr <- suppressMessages(methRead(f1, sample.id = "a", assembly = "x", pipeline = "bismarkCoverage", mincov = 1, header = FALSE))
cat("rows, default call (header = TRUE):", nrow(o_def), "\n")
cat("rows, header = FALSE:              ", nrow(o_nohdr), "\n")
cat("missing row is line 1 of the file: ", identical(getData(o_nohdr)[-1, "start"], getData(o_def)[, "start"]), "\n")
cat("first start default:", getData(o_def)$start[1], " first start header=FALSE:", getData(o_nohdr)$start[1], " file line1 start:", d1$start[1], "\n")
# the list interface used by the legacy scripts
ol <- suppressMessages(methRead(list(f1, f2), sample.id = list("a", "b"), assembly = "x", treatment = c(0, 1),
                                pipeline = "bismarkCoverage", mincov = 2))
cat("list call, mincov = 2: rows", sapply(ol, nrow), " expected with header=FALSE:", sum(d1$nC + d1$nT >= 2), sum(d2$nC + d2$nT >= 2),
    " expected if line 1 lost:", sum((d1$nC + d1$nT)[-1] >= 2), sum((d2$nC + d2$nT)[-1] >= 2), "\n")

# ---- (I-05) unknown argument high.perc ----
cat("\nformalArgs(filterByCoverage):", paste(formalArgs(filterByCoverage), collapse = ", "), "\n")
leg <- tryCatch(withCallingHandlers(
  filterByCoverage(ol, lo.count = 10, lo.perc = NULL, high.count = NULL, high.perc = 99.9),
  warning = function(w) { cat("WARNING raised:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") }),
  error = function(e) { cat("ERROR raised:", conditionMessage(e), "\n"); NULL })
good <- filterByCoverage(ol, lo.count = 10, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)
only_lo <- filterByCoverage(ol, lo.count = 10)
cat("legacy call (high.perc): rows", sapply(leg, nrow), " max coverage", sapply(leg, function(x) max(x$coverage)), "\n")
cat("lo.count only:           rows", sapply(only_lo, nrow), " max coverage", sapply(only_lo, function(x) max(x$coverage)), "\n")
cat("correct call (hi.perc):  rows", sapply(good, nrow), " max coverage", sapply(good, function(x) max(x$coverage)), "\n")
cat("legacy call identical to lo.count only:", identical(lapply(leg, getData), lapply(only_lo, getData)), "\n")
# single methylRaw object too
leg1 <- filterByCoverage(ol[[1]], lo.count = 10, lo.perc = NULL, high.count = NULL, high.perc = 99.9)
cat("methylRaw legacy call rows:", nrow(leg1), " vs lo.count only:", nrow(filterByCoverage(ol[[1]], lo.count = 10)), "\n")
# does R partial matching rescue high.perc? (supplied name must be a prefix of the formal)
f <- function(methylObj, lo.count = NULL, lo.perc = NULL, hi.count = NULL, hi.perc = NULL, ...) { list(hi.perc = hi.perc, dots = names(list(...))) }
str(f(1, lo.count = 10, lo.perc = NULL, high.count = NULL, high.perc = 99.9))

# ---- (3) operators and quantile base: methylKit against the pipeline rule ----
x <- getData(ol[[1]])                                   # as imported, mincov = 2, line 1 lost
cut_mk <- quantile(x$coverage, 99.9 / 100)              # what filterByCoverage computes (type 7, on ALL imported rows)
cat("\nmethylKit cutoff (99.9th pct of imported rows, cov >= 2):", cut_mk, " names:", names(cut_mk), "\n")
cat("99.9th pct of rows with cov >= 10 (if it were computed after lo.count):", quantile(x$coverage[x$coverage >= 10], 0.999), "\n")
manual <- x[x$coverage >= 10 & x$coverage < cut_mk, ]
cat("methylKit hi.perc result equals rows with cov >= 10 & cov < cutoff:", isTRUE(all.equal(getData(good[[1]]), manual, check.attributes = FALSE)), "\n")
cat("any kept row with coverage == floor/ceiling of cutoff? max kept:", max(good[[1]]$coverage), " cutoff:", cut_mk, "\n")
# integer cutoff case: is a CpG exactly AT the cutoff removed?
y <- ol[[1]]; yd <- getData(y)
q <- 37L
obj <- new("methylRaw", yd, sample.id = "a", assembly = "x", context = "CpG", resolution = "base")
fx <- filterByCoverage(obj, lo.count = 10, hi.count = q)
cat("hi.count = 37: max kept coverage", max(fx$coverage), "; rows with coverage == 37 in input:", sum(yd$coverage == 37), " kept:", sum(fx$coverage == 37), "\n")
cat("99.9/100 == 0.999:", 99.9/100 == 0.999, "\n")
