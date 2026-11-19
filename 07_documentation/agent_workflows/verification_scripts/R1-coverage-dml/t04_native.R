# Independent (methylKit-native) recomputation of step 04's quantities on the regional extracts.
suppressPackageStartupMessages({library(data.table); library(methylKit)})
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/R1-coverage-dml"
heads_dir <- file.path(scratch, "t05/heads"); tab <- function(f) file.path(scratch, "t04/out/04_cpg_coverage/tables", f)
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv")); setorder(lib, legacy_sample_index)
files <- as.list(file.path(heads_dir, lib$cov_file)); ids <- as.list(as.character(1:22))
sets <- list(all = 1:22, diploid = which(lib$ploidy == "diploid"), triploid = which(lib$ploidy == "triploid"), control = which(lib$treatment == "control"), heat = which(lib$treatment == "heat"))
ok <- function(label, cond) cat(sprintf("%-100s %s\n", label, if (isTRUE(cond)) "OK" else "*** FAIL ***"))
cc <- fread(tab("cpg_counts_by_coverage.csv"))
raw <- list(legacy = suppressMessages(methRead(files, sample.id = ids, assembly = "Roslin", treatment = as.integer(lib$ploidy == "triploid"), pipeline = "bismarkCoverage", mincov = 1)),
            corrected = suppressMessages(methRead(files, sample.id = ids, assembly = "Roslin", treatment = as.integer(lib$ploidy == "triploid"), pipeline = "bismarkCoverage", mincov = 1, header = FALSE)))
cut2 <- sapply(raw$corrected, function(x) quantile(x$coverage[x$coverage >= 2], 0.999))     # step 04's definition (base: coverage >= 2)
united <- list()
for (md in c("legacy", "corrected")) for (x in c(1, 3, 5, 10, 20)) {
  flt <- if (md == "legacy") filterByCoverage(raw[[md]], lo.count = x, lo.perc = NULL, high.count = NULL, high.perc = 99.9)
         else new("methylRawList", lapply(1:22, function(j) filterByCoverage(raw[[md]][[j]], lo.count = x, hi.count = unname(cut2[j]))), treatment = raw[[md]]@treatment)
  nrm <- normalizeCoverage(flt)
  for (s in names(sets)) {
    sub <- new("methylRawList", nrm@.Data[sets[[s]]], treatment = nrm@treatment[sets[[s]]])
    u <- suppressMessages(unite(sub, destrand = FALSE))
    any_n <- length(unique(unlist(lapply(sub, function(o) paste(o$chr, o$start)))))
    row <- cc[mode == md & sample_set == s & min_coverage == x]
    if (row$cpg_in_all_libraries != nrow(u) || row$cpg_in_any_library != any_n)
      cat(sprintf("MISMATCH %s %s %dx: all %d vs native %d ; any %d vs native %d\n", md, s, x, row$cpg_in_all_libraries, nrow(u), row$cpg_in_any_library, any_n))
    if (s == "all" && x %in% c(1, 10)) united[[paste(md, x)]] <- u
  }
  # methylation distribution on raw (unnormalized) counts, as step 04 does
  if (x %in% c(1, 10)) {
    dl <- fread(tab("methylation_distribution_by_library.csv"))[mode == md & min_coverage == x]
    for (j in 1:22) {
      o <- getData(flt[[j]]); h <- hist(100 * o$numCs / o$coverage, breaks = seq(0, 100, 10), plot = FALSE)$counts
      got <- dl[library_id == lib$library_id[j]][order(match(bin, c("0-10","11-20","21-30","31-40","41-50","51-60","61-70","71-80","81-90","91-100")))]$n_cpg
      if (!identical(as.integer(h), as.integer(got))) cat("MISMATCH distribution", md, x, lib$library_id[j], "\n")
    }
  }
}
cat("count and distribution comparison finished (mismatches, if any, are listed above)\n")
pv <- fread(tab("pca_variance.csv"))
for (md in c("legacy", "corrected")) for (x in c(1, 10)) {
  u <- united[[paste(md, x)]]
  d <- getData(u); mm <- as.matrix(d[, u@numCs.index]) / (as.matrix(d[, u@numCs.index]) + as.matrix(d[, u@numTs.index]))
  r_native <- cor(mm)
  r_step <- as.matrix(fread(tab(sprintf("library_correlation_%dx_%s.csv", x, md)))[, -1])
  cat(sprintf("[%s %dx] n=%d ; NaN rows native: %d ; max |r_step - r_native| = %.3g\n", md, x, nrow(d), sum(!complete.cases(mm)), max(abs(r_step - r_native))))
  pc <- PCASamples(u, obj.return = TRUE, screeplot = FALSE); dev.off()
  ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  got <- pv[mode == md & min_coverage == x][order(pc)]$percent_variance
  sc <- fread(tab(sprintf("pca_scores_%dx_%s.csv", x, md)))
  cat(sprintf("          PCA variance: max |step - native| = %.3g (PC1 native %.3f, step %.3f) ; max | |PC1 score| - |native| | = %.3g ; PC2: %.3g\n", max(abs(ve - got)), ve[1], got[1],
              max(abs(abs(sc$PC1) - abs(pc$x[, 1]))), max(abs(abs(sc$PC2) - abs(pc$x[, 2])))))
  ok(sprintf("[%s %dx] score rows labelled with the right library", md, x), identical(sc$library_id, lib$library_id) && identical(sc$ploidy, lib$ploidy) && identical(sc$treatment, lib$treatment))
}
