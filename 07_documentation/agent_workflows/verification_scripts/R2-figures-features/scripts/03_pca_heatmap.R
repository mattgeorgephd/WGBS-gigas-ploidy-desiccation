suppressPackageStartupMessages({library(data.table); library(methylKit)})
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
t05 <- file.path(repo, "04_output/05_differential_methylation/tables"); t06 <- file.path(repo, "04_output/06_dml_figures/tables")
lib <- read.csv(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv"), stringsAsFactors = FALSE); lib <- lib[order(lib$legacy_sample_index), ]
print(lib[, intersect(c("legacy_sample_index", "library_id", "ploidy", "treatment", "cohort", "acute_heat_shock", "cov_file"), names(lib))])
for (mode in c("legacy", "corrected")) {
  u <- as.data.frame(fread(file.path(t05, sprintf("united_cpg_%s.csv.gz", mode))))
  ci <- grep("^coverage", names(u)); ni <- grep("^numCs", names(u)); ti <- grep("^numTs", names(u))
  stopifnot(length(ci) == 22, all(names(u)[ci] == paste0("coverage", 1:22)))
  cat(sprintf("\n[%s] coverage == numCs + numTs in all cells: %s (cells where it fails: %d)\n", mode, all(as.matrix(u[, ci]) == as.matrix(u[, ni]) + as.matrix(u[, ti])), sum(as.matrix(u[, ci]) != as.matrix(u[, ni]) + as.matrix(u[, ti]))))
  pm <- 100 * as.matrix(u[, ni]) / (as.matrix(u[, ni]) + as.matrix(u[, ti])); colnames(pm) <- lib$library_id
  pmcov <- 100 * as.matrix(u[, ni]) / as.matrix(u[, ci])
  cat(sprintf("[%s] max |numCs/(numCs+numTs) - numCs/coverage| in percent: %.4f\n", mode, max(abs(pm - pmcov))))
  sc <- fread(file.path(t06, sprintf("pca_dml_%s.csv", mode))); va <- fread(file.path(t06, sprintf("pca_dml_variance_%s.csv", mode)))
  for (ct in c("ploidy", "heat")) {
    d <- fread(file.path(t05, sprintf("dml_%s_%s.csv", ct, mode)))
    sub <- u[d$united_row, ]
    mb <- new("methylBase", sub, sample.ids = lib$library_id, assembly = "Roslin_v1", context = "CpG", treatment = as.integer(lib[[if (ct == "ploidy") "ploidy" else "treatment"]] %in% c("triploid", "heat")),
              coverage.index = ci, numCs.index = ni, numTs.index = ti, destranded = FALSE, resolution = "base")
    pdf(NULL); pc <- PCASamples(mb, obj.return = TRUE); invisible(dev.off())
    ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
    mine <- sc[loci == paste(ct, "DML")]
    stopifnot(identical(mine$library_id, rownames(pc$x)))
    # sign of a component is arbitrary: compare up to sign
    d1 <- min(max(abs(mine$PC1 - pc$x[, 1])), max(abs(mine$PC1 + pc$x[, 1]))); d2 <- min(max(abs(mine$PC2 - pc$x[, 2])), max(abs(mine$PC2 + pc$x[, 2])))
    cat(sprintf("[%s] %s DML: methylKit::PCASamples PC1 %.4f%% PC2 %.4f%% | pipeline %.4f%% %.4f%% | max |score diff| PC1 %.3g PC2 %.3g | loci in %d, in PCA %d (pipeline says %d)\n", mode, ct, ve[1], ve[2],
                va[loci == paste(ct, "DML") & pc == 1, percent_variance], va[loci == paste(ct, "DML") & pc == 2, percent_variance], d1, d2, nrow(sub), nrow(pc$rotation), va[loci == paste(ct, "DML") & pc == 1, n_loci_used]))
    # legacy matching by start coordinate alone: how many DML starts occur on more than one sequence among the analysed CpGs?
    amb <- sum(d$start %in% u$start[duplicated(u$start)])
    first_match <- match(d$start, u$start)
    cat(sprintf("       DML whose start coordinate occurs at more than one analysed CpG: %d; DML for which the first start-only match is another row: %d\n", amb, sum(first_match != d$united_row)))
    if (sum(first_match != d$united_row) > 0) {
      sub2 <- u[first_match, ]
      mb2 <- new("methylBase", sub2, sample.ids = lib$library_id, assembly = "Roslin_v1", context = "CpG", treatment = mb@treatment, coverage.index = ci, numCs.index = ni, numTs.index = ti, destranded = FALSE, resolution = "base")
      pdf(NULL); pc2 <- PCASamples(mb2, obj.return = TRUE); invisible(dev.off()); ve2 <- 100 * pc2$sdev^2 / sum(pc2$sdev^2)
      cat(sprintf("       PCA with the legacy start-only matching (first match): PC1 %.2f%% PC2 %.2f%%\n", ve2[1], ve2[2]))
    }
    m <- pm[d$united_row, ]
    cat(sprintf("       heat map matrix %d x %d: min %.3f max %.3f; rows with zero variance: %d; any non-finite: %s\n", nrow(m), ncol(m), min(m), max(m), sum(apply(m, 1, sd) == 0), any(!is.finite(m))))
  }
  s <- fread(file.path(t05, sprintf("dml_shared_%s.csv", mode))); m <- pm[s$united_row, ]
  cat(sprintf("[%s] shared heat map matrix %d x %d: min %.3f max %.3f\n", mode, nrow(m), ncol(m), min(m), max(m)))
  cat(sprintf("[%s] shared classes: %s | sign check failures: %d | meth.diff exactly 0: %d\n", mode, paste(names(table(s$class)), table(s$class), collapse = ", "),
              sum(s$class != ifelse(s$meth.diff_ploidy > 0 & s$meth.diff_heat > 0, "both_hyper", ifelse(s$meth.diff_ploidy < 0 & s$meth.diff_heat < 0, "both_hypo", "opposite"))), sum(s$meth.diff_ploidy == 0 | s$meth.diff_heat == 0)))
  cat(sprintf("[%s] opposite: heat>0 %d, heat<0 %d; ranges ploidy diff [%.1f, %.1f], heat diff [%.1f, %.1f]\n", mode, s[class == "opposite" & meth.diff_heat > 0, .N], s[class == "opposite" & meth.diff_heat < 0, .N],
              min(s$meth.diff_ploidy), max(s$meth.diff_ploidy), min(s$meth.diff_heat), max(s$meth.diff_heat)))
  rm(u, pm, pmcov); gc(verbose = FALSE)
}
