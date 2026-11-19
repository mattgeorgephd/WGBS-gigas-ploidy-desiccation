suppressPackageStartupMessages(library(data.table))
scr <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C5-feature-odds"
cpg <- readRDS(file.path(scr, "cpg_features.rds"))
or2 <- function(dml, bg, feat) {            # step 07's table: all DML of the list vs other CpGs of the background
  a <- sum(dml & feat); b <- sum(dml & !feat); cc <- sum(bg & !dml & feat); d <- sum(bg & !dml & !feat)
  ft <- fisher.test(matrix(c(a, cc, b, d), 2)); c(ft$estimate, ft$conf.int)
}
rng <- list(c(0, 100), c(5, 95), c(10, 90), c(15, 85))
for (s in c("ploidy", "heat")) {
  dml <- cpg[[paste0("dml_", s)]]
  for (f in c("gene", "intron", "CDS", "TE", "intergenic")) {
    ft <- cpg[[paste0("in_", f)]]
    out <- sapply(rng, function(r) { x <- or2(dml, cpg$mean_ct >= r[1] & cpg$mean_ct <= r[2], ft); sprintf("%.2f (%.2f-%.2f)", x[1], x[2], x[3]) })
    cat(sprintf("%-6s %-10s | all %s | 5-95 %s | 10-90 %s | 15-85 %s\n", s, f, out[1], out[2], out[3], out[4]))
  }
  g <- sapply(rng, function(r) or2(dml, cpg$mean_ct >= r[1] & cpg$mean_ct <= r[2], cpg$in_gene)[1])
  cat(sprintf("%-6s gene-body OR as a fraction of the all-tested OR: 5-95 %.2f, 10-90 %.2f, 15-85 %.2f ; on the log scale %.2f, %.2f, %.2f\n", s,
              g[2] / g[1], g[3] / g[1], g[4] / g[1], log(g[2]) / log(g[1]), log(g[3]) / log(g[1]), log(g[4]) / log(g[1])))
}
for (r in rng) cat(sprintf("range %g-%g: n CpGs %d ; ploidy DML outside %d ; heat DML outside %d\n", r[1], r[2],
   sum(cpg$mean_ct >= r[1] & cpg$mean_ct <= r[2]), sum(cpg$dml_ploidy & !(cpg$mean_ct >= r[1] & cpg$mean_ct <= r[2])),
   sum(cpg$dml_heat & !(cpg$mean_ct >= r[1] & cpg$mean_ct <= r[2]))))
cat("tested CpGs with mean 5-10%:", sum(cpg$mean_ct >= 5 & cpg$mean_ct < 10), " of which DML (either list):", sum((cpg$dml_ploidy | cpg$dml_heat) & cpg$mean_ct >= 5 & cpg$mean_ct < 10),
    "; mean 90-95%:", sum(cpg$mean_ct > 90 & cpg$mean_ct <= 95), " DML:", sum((cpg$dml_ploidy | cpg$dml_heat) & cpg$mean_ct > 90 & cpg$mean_ct <= 95), "\n")
