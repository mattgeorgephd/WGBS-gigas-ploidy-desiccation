#!/usr/bin/env python3
"""Probe a large sorted text file on an HTTP server with range requests (no full download)."""
import subprocess, sys

def get(url, offset, length):
    end = offset + length - 1
    r = subprocess.run(["curl", "-sS", "--max-time", "120", "-r", "%d-%d" % (offset, end), url],
                       capture_output=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr.decode(errors="replace"))
    return r.stdout

def first_full_line(url, offset, size, window=4096):
    """Return (absolute offset of line start, line bytes) for the first complete line starting at or after offset."""
    if offset == 0:
        buf = get(url, 0, window)
        nl = buf.find(b"\n")
        return 0, buf[:nl]
    buf = get(url, offset - 1, min(window, size - offset + 1))
    # the byte before offset tells whether offset is a line start
    if buf[:1] == b"\n":
        start = 1
    else:
        nl = buf.find(b"\n", 1)
        if nl < 0:
            return None, None
        start = nl + 1
    nl2 = buf.find(b"\n", start)
    line = buf[start:nl2] if nl2 >= 0 else buf[start:]
    return offset - 1 + start, line

if __name__ == "__main__":
    url = sys.argv[1]; size = int(sys.argv[2]); mode = sys.argv[3]
    if mode == "sample":
        n = int(sys.argv[4])
        prev = None
        for i in range(n + 1):
            off = min(size - 2048, int(size * i / n))
            o, line = first_full_line(url, off, size)
            key = line.split(b"\t")[0].decode(errors="replace")
            flag = "" if prev is None or key >= prev else "   <-- DECREASE"
            print("%5.1f%%  offset %14d  key %s%s" % (100.0 * i / n, o, key, flag))
            prev = key
    elif mode == "search":
        target = sys.argv[4].encode()        # find smallest line-start offset whose key > target
        lo, hi = 0, size - 1
        steps = 0
        while hi - lo > 2048:
            mid = (lo + hi) // 2
            o, line = first_full_line(url, mid, size)
            key = line.split(b"\t")[0]
            steps += 1
            if key > target:
                hi = mid
            else:
                lo = mid
        print("steps", steps, "lo", lo, "hi", hi)
    elif mode == "dump":
        off = int(sys.argv[4]); length = int(sys.argv[5])
        sys.stdout.buffer.write(get(url, off, length))
