# TE track: effect of reading the RepeatMasker coordinates as BED (pipeline, legacy bedtools) against as 1-based closed
# (what the file holds), and of restricting to transposable-element classes.
suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"; t05 <- file.path(repo, "04_output/05_differential_methylation/tables"); t07 <- file.path(repo, "04_output/07_genome_features/tables")
rm <- fread(cmd = "zcat /tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/R2-figures-features/dl_ncbi_rm/rm.out.gz | tail -n +4 | awk 'BEGIN{OFS=\"\\t\"} {print $5,$6,$7,$11}'",
            header = FALSE, sep = "\t", col.names = c("chr", "b", "e", "class"))
rm[, top := sub("/.*", "", class)]
cat("records:", nrow(rm), "\n")
merged <- function(d) { d <- d[e > s]; setorder(d, chr, s, e); d[, grp := cumsum(c(TRUE, s[-1] > cummax(e)[-.N])), by = chr]; d[, .(s = min(s), e = max(e)), by = .(chr, grp)][, grp := NULL][] }
ov_one <- function(cpg, t) { out <- logical(nrow(cpg)); idx <- split(seq_len(nrow(cpg)), cpg$chr); tt <- split(t, t$chr)
  for (ch in intersect(names(idx), names(tt))) { i <- idx[[ch]]; x <- tt[[ch]]; c0 <- cpg$start[i]; k <- findInterval(c0 + 1L, x$s); hit <- k > 0; hit[hit] <- x$e[k[hit]] > c0[hit]; out[i] <- hit }
  out }
te_classes <- c("DNA", "RC", "LTR", "LINE", "SINE")
tr <- list(as_bed       = merged(rm[, .(chr, s = b, e = e)]),               # pipeline and legacy: file read as 0-based half-open
           as_1based    = merged(rm[, .(chr, s = b - 1L, e = e)]),          # what the coordinates are
           te_only      = merged(rm[top %in% te_classes, .(chr, s = b - 1L, e = e)]),
           te_or_unknown = merged(rm[top %in% c(te_classes, "Unknown"), .(chr, s = b - 1L, e = e)]),
           simple_low_sat = merged(rm[top %in% c("Simple_repeat", "Low_complexity", "Satellite", "rRNA"), .(chr, s = b - 1L, e = e)]))
for (mode in c("legacy", "corrected")) {
  u <- fread(file.path(t05, sprintf("united_cpg_%s.csv.gz", mode)))
  nl <- sum(grepl("^coverage", names(u))); mp <- rowMeans(sapply(seq_len(nl), function(j) 100 * u[[paste0("numCs", j)]] / (u[[paste0("numCs", j)]] + u[[paste0("numTs", j)]])))
  cpg <- u[, .(chr, start)]; elig <- mp >= 10 & mp <= 90; rm(u)
  dp <- fread(file.path(t05, sprintf("dml_ploidy_%s.csv", mode))); dh <- fread(file.path(t05, sprintf("dml_heat_%s.csv", mode))); ds <- fread(file.path(t05, sprintf("dml_shared_%s.csv", mode)))
  rows <- seq_len(nrow(cpg))
  sets <- list(`all analysed CpGs` = rep(TRUE, nrow(cpg)), `eligible CpGs` = elig, `ploidy DML` = rows %in% dp$united_row, `heat DML` = rows %in% dh$united_row, shared = rows %in% ds$united_row)
  ov <- sapply(tr, function(t) ov_one(cpg, t))
  pipe <- fread(file.path(t07, sprintf("feature_overlap_by_track_%s.csv", mode)))[track == "TE"]
  out <- rbindlist(lapply(names(sets), function(s) data.table(mode = mode, set = s, n = sum(sets[[s]]), pipeline_csv = pipe[set == s, n_overlap], t(colSums(ov[sets[[s]], , drop = FALSE])))))
  out[, `:=`(gained_with_1based = as_1based - as_bed, pct_as_bed = round(100 * as_bed / n, 2), pct_1based = round(100 * as_1based / n, 2), pct_te_only = round(100 * te_only / n, 2), pct_te_or_unknown = round(100 * te_or_unknown / n, 2))]
  options(width = 250); print(as.data.frame(out))
  cat(sprintf("[%s] CpGs counted as TE by the pipeline that overlap no DNA/RC/LTR/LINE/SINE record: %d of %d (%.1f%%); of those, in Unknown: %d\n", mode, sum(ov[, "as_bed"] & !ov[, "te_only"]), sum(ov[, "as_bed"]),
              100 * sum(ov[, "as_bed"] & !ov[, "te_only"]) / sum(ov[, "as_bed"]), sum(ov[, "as_bed"] & !ov[, "te_only"] & ov[, "te_or_unknown"])))
  for (b in c("all analysed CpGs", "eligible CpGs")) for (s in c("ploidy DML", "heat DML", "shared")) {
    d <- sets[[s]]; comp <- sets[[b]] & !d
    f <- function(o) { ft <- fisher.test(matrix(c(sum(d & o), sum(comp & o), sum(d & !o), sum(comp & !o)), 2)); sprintf("OR %.3f (p %.2g)", ft$estimate, ft$p.value) }
    cat(sprintf("[%s] %-18s %-10s  as BED: %s | 1-based: %s | TE classes only: %s | TE or Unknown: %s\n", mode, b, s, f(ov[, "as_bed"]), f(ov[, "as_1based"]), f(ov[, "te_only"]), f(ov[, "te_or_unknown"])))
  }
}
