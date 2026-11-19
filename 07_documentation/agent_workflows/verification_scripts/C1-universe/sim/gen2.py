# distinct lines of varying length (95 to 164 bytes), already in key order
import sys
n = int(sys.argv[1])
for i in range(1, n + 1):
    key = "XM_%09d.1" % i
    body = "%s\t105300000\tNC_047559.1\t%08d\t%08d\tNC_047559.1\tGO:%07d\t" % (key, i, i + 2, i % 9999)
    sys.stdout.write(body + "x" * (20 + (i * 37) % 70) + "\tP\n")
