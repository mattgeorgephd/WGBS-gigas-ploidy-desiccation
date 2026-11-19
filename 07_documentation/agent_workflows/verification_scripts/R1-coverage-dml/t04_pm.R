suppressPackageStartupMessages({library(data.table); library(methylKit)})
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/R1-coverage-dml"
load(file.path(scratch, "t04/out/t04_objects.RData"))
source(file.path(repo, "02_R_functions/coverage_io.R"))
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv")); setorder(lib, legacy_sample_index); n_lib <- 22; n_cpg <- nrow(cov_mat)
# the two functions exactly as in 01_pipeline/04_cpg_coverage.R
ex <- parse(file.path(repo, "01_pipeline/04_cpg_coverage.R"))
for (e in ex) if (is.call(e) && as.character(e[[1]]) == "<-" && as.character(e[[2]]) %in% c("covered", "percent_methylation")) eval(e)
files <- as.list(file.path(scratch, "t05/heads", lib$cov_file))
for (md in c("legacy", "corrected")) for (x in c(1, 10)) {
  in_all <- rep(TRUE, n_cpg); for (j in 1:22) in_all <- in_all & covered(j, x, md)
  rows <- which(in_all); pm <- percent_methylation(rows, x, md)
  pos <- key_to_position(universe[rows], chr_levels)
  raw <- suppressMessages(methRead(files, sample.id = as.list(as.character(1:22)), assembly = "R", treatment = rep(0:1, 11), pipeline = "bismarkCoverage", mincov = 1, header = md == "legacy"))
  flt <- new("methylRawList", lapply(1:22, function(j) if (md == "legacy") filterByCoverage(raw[[j]], lo.count = x) else filterByCoverage(raw[[j]], lo.count = x, hi.count = cutoffs[j])), treatment = raw@treatment)
  u <- getData(suppressMessages(unite(normalizeCoverage(flt), destrand = FALSE)))
  m <- match(paste(u$chr, u$start), paste(pos$chr, pos$start)); stopifnot(!anyNA(m), nrow(u) == length(rows))
  nat <- 100 * as.matrix(u[, seq(6, ncol(u), 3)]) / (as.matrix(u[, seq(6, ncol(u), 3)]) + as.matrix(u[, seq(7, ncol(u), 3)]))
  d <- abs(pm[m, ] - nat)
  cat(sprintf("[%s %dx] rows %d ; entries differing by > 1e-9: %d of %d (%.4f%%), in libraries: %s ; max diff %.3f points\n", md, x, length(rows), sum(d > 1e-9), length(d), 100 * mean(d > 1e-9),
              paste(which(colSums(d > 1e-9) > 0), collapse = ","), max(d)))
  w <- which(d > 1e-9, arr.ind = TRUE)
  if (nrow(w)) {
    k <- head(w, 6); for (r in seq_len(nrow(k))) {
      j <- k[r, 2]; rr <- rows[m[k[r, 1]]]; cv <- cov_mat[rr, j]; nc <- nC_mat[rr, j]
      med <- vapply(1:22, function(jj) median(cov_mat[covered(jj, x, md), jj]), numeric(1)); sf <- max(med) / med[j]; ncov <- round(sf * cv)
      cat(sprintf("    lib %d raw cov %d nC %d sf %.6f new_cov %d : step numCs=round(%s)=%d numTs=round(%s)=%d ; methylKit numCs=round(%s)=%d numTs=round(%s)=%d ; native united numCs %d numTs %d\n", j, cv, nc, sf, ncov,
                  format(ncov * nc / cv, digits = 17), round(ncov * nc / cv), format(ncov * (cv - nc) / cv, digits = 17), round(ncov * (cv - nc) / cv),
                  format(ncov * (nc / cv), digits = 17), round(ncov * (nc / cv)), format(ncov * ((cv - nc) / cv), digits = 17), round(ncov * ((cv - nc) / cv)),
                  u[k[r, 1], 6 + 3 * (j - 1)], u[k[r, 1], 7 + 3 * (j - 1)]))
    }
  }
}
