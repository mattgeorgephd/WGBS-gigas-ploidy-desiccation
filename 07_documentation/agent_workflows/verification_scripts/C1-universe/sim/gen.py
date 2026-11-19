# 24,000 distinct lines of 128 bytes, already in key order (plus a few adjacent duplicates for uniq to drop)
import sys
n = int(sys.argv[1])
for i in range(1, n + 1):
    key = "XM_%09d.1" % i
    body = "%s\t105300000\tNC_047559.1\t%08d\t%08d\tNC_047559.1\tGO:%07d\t" % (key, i, i + 2, i % 9999)
    line = body + "x" * (127 - len(body)) + "\n"
    sys.stdout.write(line)
    if i % 1000 == 0: sys.stdout.write(line)
