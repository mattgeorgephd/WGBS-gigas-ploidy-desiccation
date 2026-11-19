import sys
b = open(sys.argv[1], 'rb').read()
base = int(sys.argv[2]); nfields = int(sys.argv[3]); maxshow = int(sys.argv[4]) if len(sys.argv) > 4 else 10**9
lines = b.split(b"\n")
off = base
prev_key = None
shown = 0
for i, ln in enumerate(lines):
    key = ln.split(b"\t")[0]
    nf = ln.count(b"\t") + 1
    mark = ""
    interior = 0 < i < len(lines) - 1
    if interior and (nf != nfields or not key.startswith(b"XM_") or len(key) != 14):
        mark = "   <== MALFORMED (fields=%d)" % nf
    if prev_key is not None and key != prev_key and interior:
        mark += "   <== KEY CHANGE"
    if mark or i < 2 or i >= len(lines) - 2 or shown < maxshow:
        print("%d (mod 4096 = %4d) len %3d: %s%s" % (off, off % 4096, len(ln), ln.decode(errors="replace")[:150], mark))
        shown += 1
    off += len(ln) + 1
    prev_key = key
