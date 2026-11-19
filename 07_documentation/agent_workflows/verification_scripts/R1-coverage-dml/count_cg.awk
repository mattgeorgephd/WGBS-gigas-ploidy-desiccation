# Independent count of CG dinucleotides per sequence (case-insensitive), including pairs split across lines.
/^>/ { if (name != "") printf "%s\t%d\t%d\n", name, bp, cg; name = substr($1, 2); bp = 0; cg = 0; last = ""; next }
{
  s = toupper($0); bp += length(s)
  if (last == "C" && substr(s, 1, 1) == "G") cg++
  n = gsub(/CG/, "CG", s); cg += n
  last = substr(s, length(s), 1)
}
END { if (name != "") printf "%s\t%d\t%d\n", name, bp, cg }
