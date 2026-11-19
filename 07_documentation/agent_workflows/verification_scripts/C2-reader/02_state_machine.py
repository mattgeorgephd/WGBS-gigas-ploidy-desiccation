#!/usr/bin/env python3
"""Independent emulation of scan()'s CSV-style quote handling for a tab-separated file.

States: N (normal), A (inside an apostrophe-quoted string), D (inside a double-quoted string).
With a non-whitespace separator, scan() recognises a quote character anywhere in a field;
inside a quoted string only the same quote character closes it (a doubled quote is a literal).
A record ends at a newline met in state N. Returns the records as (first_line, last_line,
state history) so that the lines absorbed can be attributed to the quote character that held
them open.
"""
import sys, collections

def parse(path, quoteset):
    data = open(path, 'rb').read().decode('ascii')
    recs = []            # (first_line, last_line, opener_char_of_each_multiline_span)
    state = 'N'
    line = 1
    rec_start = 1
    spans = []           # (open_line, close_line, quotechar) for quotes that cross a newline
    open_line = None
    i = 0
    n = len(data)
    crossing = []        # spans of the current record that crossed a newline
    while i < n:
        c = data[i]
        if state == 'N':
            if c in quoteset:
                state = c
                open_line = line
            elif c == '\n':
                recs.append((rec_start, line, tuple(crossing)))
                crossing = []
                line += 1
                rec_start = line
        else:
            if c == state:
                # doubled quote = literal quote, stay inside
                if i + 1 < n and data[i + 1] == state:
                    i += 1
                else:
                    if open_line != line:
                        spans.append((open_line, line, state))
                        crossing.append((open_line, line, state))
                    state = 'N'
            elif c == '\n':
                line += 1
        i += 1
    eof_open = None
    if state != 'N':
        eof_open = (open_line, line - 1, state)
        spans.append((open_line, line - 1, state + ' (EOF)'))
        crossing.append((open_line, line - 1, state + ' (EOF)'))
        recs.append((rec_start, line - 1, tuple(crossing)))
    return recs, spans, eof_open, line - 1

def main():
    path = sys.argv[1]
    label = sys.argv[2]
    for qname, qs in (("defaults \"'", "\"'"), ("apostrophe only", "'"), ("double quote only", '"'), ("none", "")):
        recs, spans, eof_open, nlines = parse(path, qs)
        multi = [r for r in recs if r[1] > r[0]]
        lost = sum(r[1] - r[0] for r in recs)
        by_char = collections.Counter()
        for (a, b, ch) in spans:
            by_char[ch] += b - a
        print(f"{label:7s} quote={qname:18s} lines={nlines:6d} records={len(recs):6d} multi-line records={len(multi):3d} "
              f"lines absorbed={lost:6d} absorbed by quote char={dict(by_char)} EOF inside quote={eof_open}")
        if qs == "\"'":
            with open(f"{sys.argv[3]}/blocks_{label}.tsv", "w") as fh:
                fh.write("record_first_line\trecord_last_line\tlines_in_record\tspans(open_line,close_line,quote)\n")
                for r in multi:
                    fh.write(f"{r[0]}\t{r[1]}\t{r[1]-r[0]+1}\t{list(r[2])}\n")
            with open(f"{sys.argv[3]}/records_{label}.tsv", "w") as fh:
                fh.write("first_line\tlast_line\n")
                for r in recs:
                    fh.write(f"{r[0]}\t{r[1]}\n")
        if qs == "'":
            with open(f"{sys.argv[3]}/blocks_aposonly_{label}.tsv", "w") as fh:
                fh.write("record_first_line\trecord_last_line\tlines_in_record\tspans\n")
                for r in multi:
                    fh.write(f"{r[0]}\t{r[1]}\t{r[1]-r[0]+1}\t{list(r[2])}\n")

main()
