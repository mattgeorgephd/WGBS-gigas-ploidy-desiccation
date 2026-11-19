suppressPackageStartupMessages(library(data.table))
scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C7-legacy-go-semantics"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
P <- readRDS(file.path(scratch, "A_port.rds")); out <- P$out
d <- file.path(repo, "08_archive/tables")
h <- fread(file.path(d, "DML-getMethylDiff-HEAT-Cov10-20.csv")); p <- fread(file.path(d, "DML-getMethylDiff-PLOIDY-Cov10-20.csv"))
for (o in c("BP", "CC", "MF")) {
  m <- as.data.table(out[[paste("COMBO", o)]]$methdiff)
  m <- merge(merge(m, h[, .(chr, start, md_heat = meth.diff)], by = c("chr", "start")), p[, .(chr, start, md_ploidy = meth.diff)], by = c("chr", "start"))
  m[, class := fifelse(md_heat > 0 & md_ploidy > 0, "both positive", fifelse(md_heat < 0 & md_ploidy < 0, "both negative", "opposite signs"))]
  cat(sprintf("shared %s: %d rows; Figure 9 'hyper' (heat > 0) %d, 'hypo' (heat < 0) %d | by the three classes of Results 3.2: both positive %d, both negative %d, opposite signs %d | sign of the list value equals sign of the ploidy contrast in %d\n",
              o, nrow(m), sum(m$meth.diff > 0), sum(m$meth.diff < 0), sum(m$class == "both positive"), sum(m$class == "both negative"), sum(m$class == "opposite signs"),
              sum(sign(m$meth.diff) == sign(m$md_ploidy))))
  if (o == "BP") print(m[, .(geneID, chr, start, list_value = round(meth.diff, 2), heat = round(md_heat, 2), ploidy = round(md_ploidy, 2), class)])
}
