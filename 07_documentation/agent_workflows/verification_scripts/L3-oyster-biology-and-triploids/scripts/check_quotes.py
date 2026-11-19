#!/usr/bin/env python3
"""Count how often each intended quote occurs in the issue log (read-only) (reader L3).
Usage: python3 -I check_quotes.py <issues.md> <quotes.txt>
quotes.txt: one quote per line (lines starting with # are ignored).
"""
import sys


def main():
    with open(sys.argv[1], encoding="utf-8") as f:
        text = f.read()
    with open(sys.argv[2], encoding="utf-8") as f:
        quotes = [q.rstrip("\n") for q in f if q.strip() and not q.startswith("#")]
    for q in quotes:
        n = text.count(q)
        line = None
        if n:
            pos = text.find(q)
            line = text.count("\n", 0, pos) + 1
        print("%d x | line %s | %s" % (n, line, q[:110]))


if __name__ == "__main__":
    main()
