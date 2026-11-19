#!/bin/bash
# Stream each union file from gannet and count rows without storing it.
base="https://gannet.fish.washington.edu/panopea/WGBS-gigas-ploidy-desiccation"
for f in bed_files/BOTH/union_1x.BOTH.bedgraph.bed bed_files/BOTH/union_1x.bedgraph.bed bed_files/RONIT/union_1x.RONIT.bedgraph.bed bed_files/YAAMINI/union_1x.YAAMINI.bedgraph.bed; do
  t0=$(date +%s)
  res=$(curl -sS --fail --max-time 3000 "$base/$f" | mawk -F'\t' '{n++; nf[NF]++; if (!($1 SUBSEP $2 in seen_dummy)) {} ; w = $3-$2; wd[w]++} NR<=2{print "   head: " $0 > "/dev/stderr"} END{printf "rows=%d", n; for (k in nf) printf " NF%d=%d", k, nf[k]; for (k in wd) if (wd[k] > 0) printf " width%d=%d", k, wd[k]; printf "\n"}' 2>&1)
  t1=$(date +%s)
  echo "$f :: $res :: $((t1-t0)) s"
done
