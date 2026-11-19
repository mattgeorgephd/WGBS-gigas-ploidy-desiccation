# Print sentences of a text file that match a regular expression, with the line number of the paragraph.
# usage: kwic.py file regex [max_sentences] [context_chars]
import sys, re
fn, pat = sys.argv[1], sys.argv[2]
mx = int(sys.argv[3]) if len(sys.argv) > 3 else 40
rx = re.compile(pat, re.I)
n = 0
with open(fn, encoding="utf-8", errors="replace") as fh:
    for i, line in enumerate(fh, 1):
        if not rx.search(line):
            continue
        sents = re.split(r"(?<=[.!?])\s+(?=[A-Z(\[])", line.strip())
        for j, s in enumerate(sents):
            if rx.search(s):
                prev = sents[j-1] if j > 0 else ""
                print("L%d: %s" % (i, s))
                n += 1
                if n >= mx:
                    sys.exit(0)
