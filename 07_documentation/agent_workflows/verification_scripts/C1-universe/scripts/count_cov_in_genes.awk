# Count covered CpGs of one Bismark coverage file inside each gene.
# File 1: gene_id, chr, start, end (1-based closed). File 2: coverage (chr, 0-based start, end, pct, meth, unmeth).
# A CpG at 0-based start s occupies 1-based s+1..s+2; overlap by at least one base, any strand (as bedtools intersect).
BEGIN { FS = "\t"; OFS = "\t"; BIN = 20000 }
NR == FNR {
  n++; gid[n] = $1; gchr[n] = $2; gs[n] = $3 + 0; ge[n] = $4 + 0; cnt[n] = 0; cnt5[n] = 0
  for (b = int(gs[n] / BIN); b <= int(ge[n] / BIN); b++) { k = $2 SUBSEP b; bins[k] = (k in bins) ? bins[k] " " n : n }
  next
}
{
  total++
  p1 = $2 + 1; p2 = $2 + 2
  for (b = int(p1 / BIN); b <= int(p2 / BIN); b++) {
    k = $1 SUBSEP b
    if (k in bins) {
      m = split(bins[k], a, " ")
      for (i = 1; i <= m; i++) { j = a[i]; if (!(j in seen) && p2 >= gs[j] && p1 <= ge[j]) { cnt[j]++; if ($5 + $6 >= 5) cnt5[j]++; seen[j] = 1 } }
    }
  }
  for (j in seen) delete seen[j]
}
END {
  for (i = 1; i <= n; i++) print gid[i], gchr[i], gs[i], ge[i], cnt[i], cnt5[i]
  print "lines_in_coverage_file", total > "/dev/stderr"
}
