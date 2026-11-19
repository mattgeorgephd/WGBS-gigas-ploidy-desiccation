# Count exact occurrences of candidate quotations in the issue log (read-only). Usage: python3 -I count_quotes.py <issues.md> <quotes.txt>
import sys
text = open(sys.argv[1], encoding="utf-8").read()
for q in open(sys.argv[2], encoding="utf-8").read().split("\n"):
    q = q.rstrip("\n")
    if not q.strip():
        continue
    n = text.count(q)
    line = None
    if n:
        pos = text.find(q)
        line = text.count("\n", 0, pos) + 1
    print(n, "| line", line, "|", q[:110])
