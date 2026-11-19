# Task (a), item 6: the counts that the text takes "from a scratch computation on the cached fits".
suppressPackageStartupMessages(library(data.table)); setDTthreads(1)
options(width = 250)
R <- "/home/user/WGBS-gigas-ploidy-desiccation"
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S1-I06-numbers-and-reading"
cache <- function(a, m) readRDS(file.path(R, "06_diagnostics/dml_sensitivity/cache", sprintf("fit_%s_%s.rds", a, m)))$fit
pl <- readRDS(file.path(S, "pos_elig_legacy.rds")); pc <- readRDS(file.path(S, "pos_elig_corrected.rds"))
key <- function(p) paste(p$chr, p$start)
kl <- key(pl$pos); kc <- key(pc$pos)
in_corr <- kl %in% kc                       # legacy CpG kept by the coverage filter
cat("legacy CpGs", length(kl), "; kept by the filter", sum(in_corr), "; removed", sum(!in_corr), "; mitochondrial among tested (legacy)",
    sum(pl$pos$chr == "NC_001276.1"), "; mitochondrial kept", sum(pc$pos$chr == "NC_001276.1"), "\n")
for (contrast in c("heat", "ploidy")) {
  a <- sprintf("A1_%s_mn_f", contrast)
  fl <- cache(a, "legacy"); fc <- cache(a, "corrected")
  sig_l <- fl$qvalue < 0.01
  cat(sprintf("\n===== %s =====\n", a))
  cat("legacy: q < 0.01:", sum(sig_l), "; of these removed by the filter:", sum(sig_l & !in_corr), "; kept:", sum(sig_l & in_corr), "\n")
  cat("  removed and mitochondrial:", sum(sig_l & !in_corr & pl$pos$chr == "NC_001276.1"),
      "; removed with |difference| > 20:", sum(sig_l & !in_corr & abs(fl$meth.diff) > 20),
      "; largest |difference| among removed:", round(max(abs(fl$meth.diff[sig_l & !in_corr])), 2), "\n")
  cat("  removed by chromosome:\n"); print(table(pl$pos$chr[sig_l & !in_corr]))
  cat("  q < 0.01 legacy, mitochondrial in all (removed or kept):", sum(sig_l & pl$pos$chr == "NC_001276.1"), "\n")
  # same CpGs in both modes: p-values identical?
  i_c <- match(kl, kc)
  both <- !is.na(i_c)
  cat("  CpGs in both modes:", sum(both), "; largest |p_legacy - p_corrected|:", max(abs(fl$pvalue[both] - fc$pvalue[i_c[both]])),
      "; identical():", identical(fl$pvalue[both], fc$pvalue[i_c[both]]), "; largest |diff_legacy - diff_corrected|:",
      max(abs(fl$meth.diff[both] - fc$meth.diff[i_c[both]])), "\n")
  cat("corrected: q < 0.01:", sum(fc$qvalue < 0.01), "; of the legacy-significant CpGs kept by the filter, still q < 0.01:",
      sum(fc$qvalue[i_c[sig_l & both]] < 0.01), "; significant in corrected but not in legacy:", sum(fc$qvalue < 0.01 & !(kc %in% kl[sig_l])), "\n")
  dml_l <- sig_l & abs(fl$meth.diff) > 20
  dml_c <- fc$qvalue < 0.01 & abs(fc$meth.diff) > 20
  cat("DML legacy:", sum(dml_l), "; DML corrected:", sum(dml_c), "; corrected DML that are legacy DML:", sum(kc[dml_c] %in% kl[dml_l]), "\n")
  drop <- dml_l & !(kl %in% kc[dml_c])
  d <- data.table(chr = pl$pos$chr[drop], start = pl$pos$start[drop], in_corrected_table = in_corr[drop], p_legacy = fl$pvalue[drop],
                  p_corrected = fc$pvalue[i_c[drop]], q_legacy = fl$qvalue[drop], q_corrected = fc$qvalue[i_c[drop]], diff = fl$meth.diff[drop])
  cat("DML that leave the list in corrected mode:", nrow(d), "\n"); print(d)
  if (nrow(d) > 0) cat(sprintf("  q legacy range %.5f to %.5f; q corrected range %.5f to %.5f; p identical: %s\n", min(d$q_legacy), max(d$q_legacy),
              min(d$q_corrected, na.rm = TRUE), max(d$q_corrected, na.rm = TRUE), identical(d$p_legacy, d$p_corrected)))
  # Benjamini-Hochberg in place of SLIM
  bh_l <- p.adjust(fl$pvalue, "BH"); bh_c <- p.adjust(fc$pvalue, "BH")
  cat("BH: q < 0.01 legacy", sum(bh_l < 0.01), ", corrected", sum(bh_c < 0.01), "; DML legacy", sum(bh_l < 0.01 & abs(fl$meth.diff) > 20),
      ", corrected", sum(bh_c < 0.01 & abs(fc$meth.diff) > 20), "\n")
  r_l <- fl$qvalue / bh_l; r_c <- fc$qvalue / bh_c
  cat(sprintf("  SLIM q over BH q (the share of true nulls that SLIM estimated): legacy %.5f to %.5f; corrected %.5f to %.5f\n",
              min(r_l[bh_l < 1]), max(r_l[bh_l < 1]), min(r_c[bh_c < 1]), max(r_c[bh_c < 1])))
  print(data.table(chr = pl$pos$chr[dml_l], start = pl$pos$start[dml_l], p = fl$pvalue[dml_l], q_slim_legacy = fl$qvalue[dml_l], q_bh_legacy = bh_l[dml_l],
                   q_slim_corrected = fc$qvalue[i_c[dml_l]], q_bh_corrected = bh_c[i_c[dml_l]], diff = fl$meth.diff[dml_l]))
  # rank of the weakest call
  cat("  p-value rank of each legacy DML in legacy and in corrected mode:\n")
  print(data.table(start = pl$pos$start[dml_l], rank_legacy = rank(fl$pvalue)[dml_l], rank_corrected = rank(fc$pvalue)[i_c[dml_l]]))
}
# The chi-squared variant and the other corrected lists: stability under BH
for (a in c("A1_ploidy_mn_chisq", "A1_heat_mn_chisq", "A4_ploidy_in_heat_mn_f", "A4_ploidy_in_control_mn_f", "A2_interaction_mn_f")) {
  f <- cache(a, "legacy"); bh <- p.adjust(f$pvalue, "BH")
  cat(sprintf("%s legacy: SLIM q<0.01 %d; BH q<0.01 %d; pi0 implied %.4f\n", a, sum(f$qvalue < 0.01), sum(bh < 0.01), median((f$qvalue / bh)[bh < 1])))
}
# Uncorrected fits: SLIM share of true nulls, and what the q criterion needs
for (a in c("A0_ploidy", "A0_heat")) for (m in c("legacy", "corrected")) {
  f <- cache(a, m); bh <- p.adjust(f$pvalue, "BH")
  cat(sprintf("%s %s: SLIM q<0.01 %d; BH q<0.01 %d; pi0 implied %.4f; largest p with q<0.01 %.3g; CpGs above 20 points without q<0.01: %d\n", a, m,
              sum(f$qvalue < 0.01), sum(bh < 0.01), median((f$qvalue / bh)[bh < 1]), max(f$pvalue[f$qvalue < 0.01]), sum(abs(f$meth.diff) > 20 & !(f$qvalue < 0.01))))
  w <- which(abs(f$meth.diff) > 20 & !(f$qvalue < 0.01))
  if (length(w)) { p <- if (m == "legacy") pl$pos else pc$pos; print(cbind(p[w], f[w])) }
}
