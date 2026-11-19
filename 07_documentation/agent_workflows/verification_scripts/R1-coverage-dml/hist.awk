BEGIN { FS = "\t"; bad_end = 0; unsorted = 0; nf_bad = 0; prev_chr = ""; prev_start = -1; nchr = 0; dupchr = 0 }
{
  if (NF != 6) nf_bad++
  c = $5 + $6
  h[c]++
  if ($3 != $2 + 2) bad_end++
  if ($1 != prev_chr) { if ($1 in seen) dupchr++; seen[$1] = 1; nchr++; prev_chr = $1; prev_start = -1 }
  if ($2 + 0 <= prev_start) unsorted++
  prev_start = $2 + 0
  if (NR == 1) first = $0
  if (c >= 10) { n10++; }
}
END {
  for (k in h) printf "hist\t%d\t%d\n", k, h[k]
  printf "nlines\t%d\n", NR
  printf "bad_end\t%d\n", bad_end
  printf "nf_bad\t%d\n", nf_bad
  printf "unsorted_within_chr\t%d\n", unsorted
  printf "chr_blocks\t%d\n", nchr
  printf "chr_reappears\t%d\n", dupchr
  printf "first\t%s\n", first
}
