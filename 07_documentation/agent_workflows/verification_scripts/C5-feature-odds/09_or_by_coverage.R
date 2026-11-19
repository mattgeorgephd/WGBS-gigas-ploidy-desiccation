suppressPackageStartupMessages(library(data.table))
scr <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C5-feature-odds"
cpg <- readRDS(file.path(scr, "cpg_features.rds"))
q <- quantile(cpg$mean_cov[cpg$elig_p], 0:4 / 4)
cat("quartile cut points of mean coverage among eligible CpGs:", round(q, 2), "\n")
or2 <- function(dml, bg, feat) {
  a <- sum(dml & feat & bg); b <- sum(dml & !feat & bg); cc <- sum(bg & !dml & feat); d <- sum(bg & !dml & !feat)
  ft <- fisher.test(matrix(c(a, cc, b, d), 2))
  sprintf("%.2f (%.2f-%.2f) [%d/%d vs %d/%d]", ft$estimate, ft$conf.int[1], ft$conf.int[2], a, a + b, cc, cc + d)
}
for (s in c("dml_ploidy", "dml_heat")) for (f in c("gene", "intron", "CDS", "TE", "intergenic")) {
  dml <- cpg[[s]]; ft <- cpg[[paste0("in_", f)]]
  cat(sprintf("%-10s %-10s all eligible %s | cov <= Q3 %s | cov > Q3 %s\n", s, f,
              or2(dml, cpg$elig_p, ft), or2(dml, cpg$elig_p & cpg$mean_cov <= q[4], ft), or2(dml, cpg$elig_p & cpg$mean_cov > q[4], ft)))
}
cat("\nshare of eligible CpGs in the top coverage quartile, by feature:\n")
e <- cpg[elig_p == TRUE]
for (f in c("gene", "intron", "CDS", "TE", "intergenic")) cat(sprintf("%-10s in: %.1f%%  out: %.1f%%\n", f,
   100 * mean(e$mean_cov[e[[paste0("in_", f)]]] > q[4]), 100 * mean(e$mean_cov[!e[[paste0("in_", f)]]] > q[4])))
