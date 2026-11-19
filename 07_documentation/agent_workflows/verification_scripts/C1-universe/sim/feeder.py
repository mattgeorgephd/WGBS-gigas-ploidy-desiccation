# Pass stdin to stdout slowly; at a given line number, pause and delete a ballast file (frees disk space).
import sys, os, time
trigger = int(sys.argv[1]); ballast = sys.argv[2] if len(sys.argv) > 2 else None
out = sys.stdout.buffer
for i, line in enumerate(sys.stdin.buffer, 1):
    out.write(line)
    if i % 200 == 0:
        out.flush(); time.sleep(0.002)
    if ballast and i == trigger:
        out.flush(); time.sleep(1.5)
        os.remove(ballast); sys.stderr.write("[feeder] ballast removed after line %d\n" % i)
        time.sleep(0.5)
out.flush()
