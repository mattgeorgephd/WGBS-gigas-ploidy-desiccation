# venkataraman_2022_check.R
#
# Purpose   The DML call of Venkataraman et al. 2022 (BMC Genomics 23:556,
#           doi:10.1186/s12864-022-08781-5), the study whose Methods the draft
#           follows (I-06, I-28), reproduced from its own data and code and set
#           against every regrouping of its libraries (author's request of
#           2026-11-26). Their code is RobertsLab/project-gigas-oa-meth (GitHub;
#           mirrored on OSF, doi:10.17605/OSF.IO/YGCTB), code/06-methylKit.Rmd at
#           commit 77463bbbebc8f4258ca6b39c8a187e9dac4be169; their coverage files are
#           on gannet. Four parts:
#             1. Import as their script does: methRead (bismarkCoverage, mincov = 2),
#                filterByCoverage(lo.count = 5, lo.perc = NULL, high.count = NULL,
#                high.perc = 99.9), normalizeCoverage, unite(destrand = FALSE). The
#                arguments high.count and high.perc are not arguments of
#                filterByCoverage (hi.count and hi.perc in methylKit 1.17.4, the
#                version of the paper, source commit 9620cd6 of al2na/methylKit, and
#                in 1.28.0): they fall into `...` and are ignored, so the 99.9th
#                percentile filter that the paper describes was not applied, and it
#                is not applied here (issues.md, I-05 and I-06).
#             2. Their model: calculateDiffMeth(covariates = gonad stage,
#                overdispersion = "MN", test = "Chisq"), DML at q < 0.01 and a
#                difference of more than 50 points (getMethylDiff). Fitted in blocks
#                of rows (resumable); the SLIM q-values are computed once over all
#                CpGs with methylKit's p.adjusted(method = "SLIM"). The DML are checked
#                against their list DML-pH-50-Cov5-All.csv (1,599).
#             3. Every split of the 8 libraries into 4 and 4 (35 up to the mirror
#                image; the real one among them), the stage of each library kept.
#                The test is methylKit 1.17.4's logReg (R/diffMeth.R) written out in
#                base R; for the real labels it is checked against part 2. Per
#                labelling: CpGs whose difference exceeds 50 points with a corrected
#                p at or below the threshold of the real list; those below 0.01 /
#                (CpGs tested); and the share of small p-values at two seeded samples
#                of 3,000 eligible CpGs (mean 25% to 75%).
#             4. Summary.
#           The stage covariate is a character column, a factor of 4 levels in the
#           model, so the full model has 5 parameters for 8 libraries: the dispersion
#           of each CpG is estimated from 3 residual degrees of freedom. One split
#           cannot be told apart from stage and is not estimable.
#           Their removal of DML at C-to-T SNPs (315 of 1,599) comes after the call and
#           is not part of either side of the comparison.
#           Not in diagnostics.csv: it needs the network and data of another study.
#           Nothing in 01_pipeline/ reads its outputs.
#
# Run       WGBS_CORES=4 Rscript 06_diagnostics/venkataraman_2022_check.R
#           (downloads 3.7 GB into the download cache, subfolder venkataraman_2022;
#           about 10 minutes for part 1 and 60 for part 2 on 4 cores, 2 for part 3;
#           peak memory about 9 GB)
#
# Writes    06_diagnostics/venkataraman_2022_check/
#             reproduction.csv   their DML against the ones reproduced here
#             relabel.csv        one row per split of the libraries
#             summary.csv        the headline numbers
#             session_info.txt
#           06_diagnostics/venkataraman_2022_check/cache/ (git-ignored): united.rds,
#             parts/part_<k>.rds, fit.rds

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "methylKit"))
t0 <- Sys.time()
log_msg("venkataraman_2022_check: the DML call of Venkataraman et al. 2022, reproduced and set against every split of its libraries")

their_commit <- "77463bbbebc8f4258ca6b39c8a187e9dac4be169"
gannet_url <- "https://gannet.fish.washington.edu/spartina/project-gigas-oa-meth/output/bismark-roslin/"
their_dml_url <- sprintf("https://raw.githubusercontent.com/RobertsLab/project-gigas-oa-meth/%s/output/06-methylKit/DML/DML-pH-50-Cov5-All.csv",
                         their_commit)
cores <- max(1L, as.integer(Sys.getenv("WGBS_CORES", "1")))
out_dir <- here("06_diagnostics", "venkataraman_2022_check")
cache <- file.path(out_dir, "cache"); parts <- file.path(cache, "parts")
dir.create(parts, recursive = TRUE, showWarnings = FALSE)
raw_dir <- gannet_cache_dir(config, "venkataraman_2022")
out_file <- function(name) file.path(out_dir, name)
fetch <- function(url, dest) {
  if (file.exists(dest) && file.size(dest) > 0) return(dest)
  log_msg("downloading %s", basename(dest))
  if (utils::download.file(url, dest, mode = "wb", quiet = TRUE) != 0) stop("Download failed: ", url, call. = FALSE)
  dest
}

# Their design (code/06-methylKit.Rmd): libraries 1 to 4 low pH (1), 5 to 8 ambient (0);
# covariateMetadataAll: gonad stage.
treatment <- c(1, 1, 1, 1, 0, 0, 0, 0)
stage_chr <- c("2", "0", "1", "1", "0", "1", "2", "3")
files <- sprintf("zr3616_%d_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov", 1:8)

# ---- 1. Import ---------------------------------------------------------------------------------
united_file <- file.path(cache, "united.rds")
if (!file.exists(united_file)) {
  paths <- vapply(files, function(f) fetch(paste0(gannet_url, f), file.path(raw_dir, f)), "")
  filtered <- lapply(seq_along(paths), function(i) {
    r <- methRead(paths[[i]], sample.id = as.character(i), assembly = "Roslin", treatment = treatment[i],
                  pipeline = "bismarkCoverage", mincov = 2)
    # Their call, verbatim; high.count and high.perc are ignored by methylKit (see the header).
    f <- filterByCoverage(r, lo.count = 5, lo.perc = NULL, high.count = NULL, high.perc = 99.9)
    log_msg("  library %d: %d CpGs at 2x, %d after the filter", i, nrow(r), nrow(f))
    rm(r); invisible(gc()); f
  })
  rl <- normalizeCoverage(new("methylRawList", filtered, treatment = treatment)); rm(filtered); invisible(gc())
  united <- unite(rl, destrand = FALSE, mc.cores = min(cores, 2L)); rm(rl); invisible(gc())
  saveRDS(united, united_file)
}
united <- readRDS(united_file)
n <- nrow(united)
log_msg("part 1: %d CpGs with 5x in all 8 libraries", n)

# ---- 2. Their model ------------------------------------------------------------------------------
cov <- data.frame(stage = stage_chr)
block <- 250000L; starts <- seq(1L, n, by = block)
for (k in seq_along(starts)) {
  f <- file.path(parts, sprintf("part_%03d.rds", k))
  if (file.exists(f)) next
  rows <- starts[k]:min(n, starts[k] + block - 1L)
  d <- suppressWarnings(calculateDiffMeth(methylKit::select(united, rows), covariates = cov, overdispersion = "MN",
                                          test = "Chisq", mc.cores = cores))
  x <- getData(d)
  saveRDS(data.table(row = rows, chr = x$chr, start = x$start, pvalue = x$pvalue, meth.diff = x$meth.diff), f)
  log_msg("  block %d of %d", k, length(starts))
}
fit <- rbindlist(lapply(file.path(parts, sprintf("part_%03d.rds", seq_along(starts))), readRDS))
stopifnot(identical(fit$row, seq_len(n)))
p_adjusted <- utils::getFromNamespace("p.adjusted", "methylKit")
fit[, qvalue := p_adjusted(pvalue, method = "SLIM")]
saveRDS(fit, file.path(cache, "fit.rds"))
dml <- fit[qvalue < 0.01 & abs(meth.diff) > 50]
theirs <- fread(fetch(their_dml_url, file.path(raw_dir, "DML-pH-50-Cov5-All.csv")))
m <- merge(fit, theirs[, .(chr, start, p_theirs = pvalue, d_theirs = meth.diff)], by = c("chr", "start"))
p_star <- max(fit$pvalue[fit$qvalue < 0.01])
repro <- data.table(cpgs_united = n, dml_here = nrow(dml), dml_theirs = nrow(theirs),
                    dml_in_both = length(intersect(paste(dml$chr, dml$start), paste(theirs$chr, theirs$start))),
                    their_dml_among_cpgs_here = nrow(m),
                    largest_relative_p_difference = max(abs(m$pvalue - m$p_theirs) / m$p_theirs),
                    largest_meth_diff_difference = max(abs(m$meth.diff - m$d_theirs)),
                    p_threshold_of_list = p_star, their_code_commit = their_commit)
write_table(repro, out_file("reproduction.csv"))
log_msg("part 2: %d DML here, %d theirs, %d in both; largest relative p difference %.1e", repro$dml_here, repro$dml_theirs,
        repro$dml_in_both, repro$largest_relative_p_difference)
if (repro$dml_in_both != repro$dml_theirs || repro$dml_here != repro$dml_theirs) stop("The DML do not reproduce their list.", call. = FALSE)

# ---- 3. Every split of the libraries -----------------------------------------------------------
u <- as.data.table(getData(united)); rm(united); invisible(gc())
cs <- as.matrix(u[, paste0("numCs", 1:8), with = FALSE]) + 0
ns <- cs + as.matrix(u[, paste0("numTs", 1:8), with = FALSE])
stopifnot(identical(u$start, fit$start)); rm(u); invisible(gc())
stage <- factor(stage_chr)
p_single <- 0.01 / n
quiet_glm <- function(X, y, w) suppressWarnings(stats::glm.fit(X, y / w, weights = w, family = stats::binomial()))
one_test <- function(y, w, Xf, Xr) {                    # logReg of methylKit 1.17.4, MN, chi-squared
  full <- quiet_glm(Xf, y, w); red <- quiet_glm(Xr, y, w)
  mu <- full$fitted.values
  phi <- sum((y - w * mu)^2 / (mu * (w - w * mu))) / (length(w) - ncol(Xf))
  dev <- red$deviance - full$deviance
  p_unc <- pchisq(dev, 1, lower.tail = FALSE); p_cor <- pchisq(dev / ifelse(phi > 1, phi, 1), 1, lower.tail = FALSE)
  c(p_unc = if (is.na(p_unc)) 1 else p_unc, p_cor = if (is.na(p_cor)) 1 else p_cor, phi = phi)
}
run <- function(rows, Xf, Xr) if (!length(rows)) matrix(numeric(0), 0, 3) else
  t(vapply(rows, function(i) one_test(cs[i, ], ns[i, ], Xf, Xr), numeric(3)))
pooled_diff <- function(z) 100 * (rowSums(cs[, z == 1]) / rowSums(ns[, z == 1]) - rowSums(cs[, z == 0]) / rowSums(ns[, z == 0]))
sets <- combn(8, 4); sets <- sets[, sets[1, ] == 1, drop = FALSE]
labels <- apply(sets, 2, function(s) { z <- integer(8); z[s] <- 1L; z })
is_real <- apply(labels, 2, function(z) identical(as.numeric(z), treatment))
stopifnot(ncol(labels) == 35, sum(is_real) == 1)
labels <- cbind(labels[, is_real], labels[, !is_real])
Xr <- model.matrix(~ stage)
mean_pm <- rowMeans(100 * cs / ns); eligible <- which(mean_pm >= 25 & mean_pm <= 75)
set.seed(606); s1 <- sort(sample(eligible, 3000)); set.seed(607); s2 <- sort(sample(setdiff(eligible, s1), 3000))
res <- parallel::mclapply(seq_len(ncol(labels)), function(j) {
  z <- labels[, j]
  Xf <- model.matrix(~ z + stage)
  estimable <- qr(Xf)$rank == ncol(Xf)
  d <- pooled_diff(z); above <- which(abs(d) > 50)
  ta <- if (estimable) run(above, Xf, Xr) else matrix(NA_real_, length(above), 3)
  lev <- if (estimable) rbind(run(s1, Xf, Xr), run(s2, Xf, Xr)) else matrix(NA_real_, 6000, 3)
  called <- ta[, 2] <= p_star * (1 + 1e-9)
  out <- data.table(labelling = j - 1L, real = j == 1L, low_ph_group = paste(which(z == 1), collapse = ","), estimable = estimable,
                    n_above_50_points = length(above),
                    n_called = if (estimable) sum(called) else NA_integer_,
                    n_called_dispersion_not_above_1 = if (estimable) sum(called & ta[, 3] <= 1, na.rm = TRUE) else NA_integer_,
                    n_called_on_their_own = if (estimable) sum(ta[, 2] < p_single) else NA_integer_,
                    share_corrected_p_below_0_05 = mean(lev[, 2] < 0.05), share_corrected_p_below_0_01 = mean(lev[, 2] < 0.01),
                    share_corrected_p_below_0_001 = mean(lev[, 2] < 0.001), share_uncorrected_p_below_0_01 = mean(lev[, 1] < 0.01),
                    median_dispersion_eligible = stats::median(lev[, 3], na.rm = TRUE),
                    share_dispersion_not_above_1_eligible = mean(lev[, 3] <= 1, na.rm = TRUE))
  if (j == 1L) attr(out, "check") <- c(p = max(abs(ta[, 2] - fit$pvalue[above]) / pmax(fit$pvalue[above], .Machine$double.xmin)),
                                       d = max(abs(d - fit$meth.diff)))
  out
}, mc.cores = cores, mc.preschedule = FALSE)
bad <- vapply(res, function(r) inherits(r, "try-error"), logical(1))
if (any(bad)) stop("Split ", which(bad)[1], " failed: ", as.character(res[[which(bad)[1]]]), call. = FALSE)
chk <- attr(res[[1]], "check"); r <- rbindlist(res)
if (chk[["p"]] > 1e-6 || chk[["d"]] > 1e-8 || r[real == TRUE, n_called] != nrow(dml)) {
  stop("Part 3 does not reproduce the real fit (p ", signif(chk[["p"]], 3), ", difference ", signif(chk[["d"]], 3), ", called ",
       r[real == TRUE, n_called], " of ", nrow(dml), ").", call. = FALSE)
}
write_table(r, out_file("relabel.csv"))

# ---- 4. Summary ----------------------------------------------------------------------------------
o <- r[real == TRUE]; x <- r[real == FALSE & estimable == TRUE]
summ <- data.table(
  measure = c("DML (q < 0.01, more than 50 points)", "DML with p below 0.01 / CpGs tested",
              "share of corrected p below 0.05 at eligible CpGs", "share of corrected p below 0.01 at eligible CpGs",
              "share of corrected p below 0.001 at eligible CpGs", "DML at which the dispersion estimate is not above 1"),
  real = c(o$n_called, o$n_called_on_their_own, o$share_corrected_p_below_0_05, o$share_corrected_p_below_0_01,
           o$share_corrected_p_below_0_001, o$n_called_dispersion_not_above_1),
  other_splits_median = c(stats::median(x$n_called), stats::median(x$n_called_on_their_own), stats::median(x$share_corrected_p_below_0_05),
                          stats::median(x$share_corrected_p_below_0_01), stats::median(x$share_corrected_p_below_0_001),
                          stats::median(x$n_called_dispersion_not_above_1)),
  other_splits_min = c(min(x$n_called), min(x$n_called_on_their_own), min(x$share_corrected_p_below_0_05),
                       min(x$share_corrected_p_below_0_01), min(x$share_corrected_p_below_0_001), min(x$n_called_dispersion_not_above_1)),
  other_splits_max = c(max(x$n_called), max(x$n_called_on_their_own), max(x$share_corrected_p_below_0_05),
                       max(x$share_corrected_p_below_0_01), max(x$share_corrected_p_below_0_001), max(x$n_called_dispersion_not_above_1)),
  other_splits_at_or_above_real = c(sum(x$n_called >= o$n_called), sum(x$n_called_on_their_own >= o$n_called_on_their_own),
                                    sum(x$share_corrected_p_below_0_05 >= o$share_corrected_p_below_0_05),
                                    sum(x$share_corrected_p_below_0_01 >= o$share_corrected_p_below_0_01),
                                    sum(x$share_corrected_p_below_0_001 >= o$share_corrected_p_below_0_001), NA),
  other_splits_estimable = nrow(x), splits_not_estimable = sum(!r$estimable))
summ[, permutation_p := fifelse(is.na(other_splits_at_or_above_real), NA_real_, (other_splits_at_or_above_real + 1) / (other_splits_estimable + 1))]
write_table(summ, out_file("summary.csv"))
print(summ)
writeLines(capture.output(sessionInfo()), out_file("session_info.txt"))
log_msg("venkataraman_2022_check finished in %.1f min; outputs in %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), out_dir)
