# tests/check_overdispersion.R
#
# T3-c of the settings layout (settings_build_spec.md, section 8): after
#   Rscript run_pipeline.R --overdispersion
# the DML lists of the analysis 'mnchisq' (04_output/not_adopted/mnchisq/) must
# be exactly the lists of 06_diagnostics/dml_sensitivity/ for the overdispersion
# correction with the chi-squared test on the paper's CpG set (analysis A1):
# 14 ploidy and 185 heat DML, at the same positions, with the same p and q
# (to 1e-9, relative), every one of them in the paper's lists; 480,579 CpGs tested.
# Run from the repository root; exit status 1 when a check fails.

suppressPackageStartupMessages(library(data.table))
na <- function(f) file.path("04_output/not_adopted/mnchisq/05_differential_methylation/tables", f)
paper <- function(f) file.path("04_output/05_differential_methylation/tables", f)
pass <- TRUE
for (contrast in c("ploidy", "heat")) {
  d <- fread(na(sprintf("dml_%s_mnchisq.csv", contrast)))
  ref <- fread(sprintf("06_diagnostics/dml_sensitivity/dml_A1_%s_mn_chisq_legacy.csv", contrast))
  p <- fread(paper(sprintf("dml_%s_legacy.csv", contrast)))
  same_pos <- identical(paste(d$chr, d$start), paste(ref$chr, ref$start))
  same_pq <- same_pos && isTRUE(all.equal(d$pvalue, ref$pvalue, tolerance = 1e-9)) && isTRUE(all.equal(d$qvalue, ref$qvalue, tolerance = 1e-9))
  in_paper <- all(paste(d$chr, d$start) %in% paste(p$chr, p$start))
  ok <- same_pos && same_pq && in_paper
  cat(sprintf("%-4s %s: %d DML (diagnostic %d); same positions %s, same p and q %s, all in the paper's list (%d) %s\n",
              if (ok) "ok" else "FAIL", contrast, nrow(d), nrow(ref), same_pos, same_pq, nrow(p), in_paper))
  pass <- pass && ok
}
h <- fread(na("dml_counts_by_mode.csv"))
ok <- h$cpg_analysed == 480579 && h$mode == "mnchisq"
cat(sprintf("%-4s CpGs tested: %d; shared DML: %d\n", if (ok) "ok" else "FAIL", h$cpg_analysed, h$shared))
quit(status = if (pass && ok) 0 else 1)
