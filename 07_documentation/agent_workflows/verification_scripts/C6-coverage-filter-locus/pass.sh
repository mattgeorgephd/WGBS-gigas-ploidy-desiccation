#!/bin/bash
# One mawk pass per coverage file, one file at a time.
# Per library: coverage histogram, rows in NC_047563.1 +-100 kb of the locus, rows at the
# first-line positions of any file, rows with coverage >= 50, and the running intersection
# of CpGs with >= 10 reads (complete read, no upper filter), carrying cov:nC per library.
S=/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C6-coverage-filter-locus
C=/home/user/gannet_cache/cov_files
order="zr3644_14 zr3644_3 zr3644_5 zr3644_16 zr3644_4 zr3644_1 zr3644_18 zr3644_6 zr3644_13 zr3644_2 zr3644_17 zr3644_15 zr3534_7 zr3534_9 zr3534_4 zr3534_6 zr3534_2 zr3534_10 zr3534_1 zr3534_3 zr3534_5 zr3534_8"
mkdir -p $S/pass
echo "$order" | tr ' ' '\n' > $S/pass/order.txt
cut -f2,3 $S/first_lines.tsv | sort -u > $S/pass/first_positions.tsv
: > $S/pass/progress.log
prev=""
for lib in $order; do
  f=$C/${lib}_R1.CpG_report.merged_CpG_evidence.cov
  nice -n 10 mawk -F'\t' -v OFS='\t' -v lib=$lib -v out=$S/pass -v prev="$prev" -v fpos=$S/pass/first_positions.tsv '
  BEGIN {
    while ((getline line < fpos) > 0) { split(line, a, "\t"); fp[a[1] SUBSEP a[2]] = 1 }
    close(fpos)
    if (prev != "") {
      while ((getline line < prev) > 0) { split(line, a, "\t"); cand[a[1] SUBSEP a[2]] = a[3] }
      close(prev)
    }
    newc = out "/cand_after_" lib ".tsv"
    reg = out "/region_" lib ".tsv"; fpo = out "/fpos_" lib ".tsv"; hi = out "/hi_" lib ".tsv"
    printf "" > newc; printf "" > reg; printf "" > fpo; printf "" > hi
  }
  {
    cov = $5 + $6; h[cov]++
    if (cov >= 10) {
      if (prev == "") print $1, $2, cov ":" $5 > newc
      else { k = $1 SUBSEP $2; if (k in cand) print $1, $2, cand[k] "," cov ":" $5 > newc }
      if (cov >= 50) print $1, $2, $5, $6 > hi
    }
    if ($2 <= 16536 && (($1 SUBSEP $2) in fp)) print lib, $0 > fpo
    if ($1 == "NC_047563.1" && $2 >= 39533323 && $2 <= 39735520) print lib, $0 > reg
  }
  END {
    hf = out "/hist_" lib ".tsv"
    for (c in h) print c, h[c] > hf
    print lib, NR > (out "/nlines_" lib ".txt")
  }' "$f"
  new=$S/pass/cand_after_${lib}.tsv
  echo "$(date +%T) $lib done: $(wc -l < $new) candidates" >> $S/pass/progress.log
  if [ -n "$prev" ]; then rm -f "$prev"; fi
  prev=$new
done
echo "ALL DONE $(date +%T)" >> $S/pass/progress.log
