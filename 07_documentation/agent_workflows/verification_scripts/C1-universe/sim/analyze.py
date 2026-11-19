import sys, re
b = open(sys.argv[1], 'rb').read()
print("file size", len(b), "bytes; = 4096 x %.3f" % (len(b) / 4096))
lines = b.split(b"\n")
off = 0; prev = None; bad = 0; jumps = []
for ln in lines[:-1]:
    m = re.match(rb"^XM_(\d{9})\.1\t", ln)
    ok = m is not None and len(ln) == 127
    if not ok:
        bad += 1
        # where does a second line start inside this malformed line?
        k = ln.find(b"XM_", 1)
        print("  malformed line at offset %d (len %d): %r" % (off, len(ln), ln[:60] + b" ... " + ln[-40:]))
        if k > 0:
            print("    a fresh line starts at offset %d = 4096 x %.3f ; text there: %r" % (off + k, (off + k) / 4096, ln[k:k+20]))
    key = int(m.group(1)) if m else None
    if key is not None and prev is not None and key != prev + 1:
        jumps.append((off, prev, key))
    if key is not None: prev = key
    off += len(ln) + 1
print("lines", len(lines) - 1, "malformed", bad, "; last key", prev, "; ends with newline:", b.endswith(b"\n"))
for o, a, c in jumps: print("  key jump at offset %d (4096 x %.3f): %d -> %d (missing %d keys)" % (o, o / 4096, a, c, c - a - 1))
