# Streaming count of CG dinucleotides in a FASTA, handling pairs that span a line break.
# Per sequence: length, CG (case-insensitive), CG with both bases upper case, both lower case, mixed case,
# number of non-ACGT characters.
function flush() {
  if (name != "") printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", name, len, cg_any, cg_up, cg_low, cg_mixed, nonacgt
}
/^>/ {
  flush()
  name = substr($1, 2); len = 0; cg_any = 0; cg_up = 0; cg_low = 0; cg_mixed = 0; nonacgt = 0; prev = ""
  next
}
{
  line = $0
  len += length(line)
  s = prev line                       # prepend last base of the previous line
  u = toupper(s)
  t = u; cg_any += gsub(/CG/, "", t)
  t = s; n_up = gsub(/CG/, "", t); cg_up += n_up
  t = s; n_low = gsub(/cg/, "", t); cg_low += n_low
  t = s; n_m1 = gsub(/Cg/, "", t); t = s; n_m2 = gsub(/cG/, "", t); cg_mixed += n_m1 + n_m2
  t = toupper(line); nonacgt += gsub(/[^ACGT]/, "", t)
  prev = substr(line, length(line), 1)
}
END { flush() }
