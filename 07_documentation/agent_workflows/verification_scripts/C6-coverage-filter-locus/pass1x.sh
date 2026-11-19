#!/bin/bash
# CpGs with at least 1 read in all 22 libraries (complete read of every file): running intersection, one file at a time.
S=/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C6-coverage-filter-locus
C=/home/user/gannet_cache/cov_files
order="zr3644_14 zr3644_1 zr3644_3 zr3644_17 zr3644_2 zr3644_16 zr3644_5 zr3644_18 zr3644_4 zr3644_6 zr3534_4 zr3534_2 zr3534_1 zr3644_13 zr3534_5 zr3644_15 zr3534_3 zr3534_7 zr3534_9 zr3534_10 zr3534_6 zr3534_8"
mkdir -p $S/pass1x; : > $S/pass1x/progress.log
prev=""
for lib in $order; do
  f=$C/${lib}_R1.CpG_report.merged_CpG_evidence.cov
  new=$S/pass1x/cand_${lib}.tsv
  if [ -z "$prev" ]; then
    nice -n 15 cut -f1,2 "$f" > "$new"
  else
    nice -n 15 mawk -F'\t' 'NR==FNR{c[$1 FS $2]; next} (($1 FS $2) in c){print $1 FS $2}' "$prev" "$f" > "$new"
    rm -f "$prev"
  fi
  echo "$(date +%T) $lib: $(wc -l < $new)" >> $S/pass1x/progress.log
  prev=$new
done
mv "$prev" $S/pass1x/final_1x_all22.tsv
echo "ALL DONE $(date +%T)" >> $S/pass1x/progress.log
