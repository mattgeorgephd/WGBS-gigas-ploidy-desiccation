import sys, re
b = open(sys.argv[1], 'rb').read()
print("file size", len(b), "bytes")
off = 0; prev = None; bad = 0
for ln in b.split(b"\n")[:-1]:
    m = re.match(rb"^XM_(\d{9})\.1\t[^\n]*\tP$", ln)
    wellformed = m is not None and ln.count(b"XM_") == 1 and ln.count(b"\t") == 8
    if not wellformed:
        bad += 1
        k = ln.find(b"XM_", 1)
        print("  malformed line at offset %d, len %d, fields %d: %r" % (off, len(ln), ln.count(b"\t") + 1, ln[:70] + b"..." ))
        if k > 0:
            m2 = re.match(rb"XM_(\d{9})", ln[k:])
            print("    cut after %d bytes; a fresh line (key %d) starts at offset %d = 4096 x %.4f" % (k, int(m2.group(1)), off + k, (off + k) / 4096))
            m1 = re.match(rb"^XM_(\d{9})", ln)
            print("    keys: %d -> %d (missing %d keys in between)" % (int(m1.group(1)), int(m2.group(1)), int(m2.group(1)) - int(m1.group(1)) - 1))
    off += len(ln) + 1
last = b.rstrip(b"\n").split(b"\n")[-1]
print("malformed lines:", bad, "; last line key:", last[:14].decode(), "; ends with newline:", b.endswith(b"\n"))
