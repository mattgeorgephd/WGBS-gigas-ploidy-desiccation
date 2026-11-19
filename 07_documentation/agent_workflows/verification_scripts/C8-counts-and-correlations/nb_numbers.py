import json, sys, re, glob
rx = re.compile(r"(?<![\d.])(1[0-4][\s,]?\d{3}[\s,]?\d{3})(?![\d.])")
for p in sorted(glob.glob("**/*.ipynb", recursive=True)):
    nb = json.load(open(p, encoding="utf-8"))
    cells = nb.get("cells") or [c for ws in nb.get("worksheets", []) for c in ws.get("cells", [])]
    print(f"=== {p}: {len(cells)} cells")
    for i, c in enumerate(cells):
        src = "".join(c.get("source", c.get("input", "")))
        outs = []
        for o in c.get("outputs", []):
            if "text" in o: outs.append("".join(o["text"]))
            for k, v in (o.get("data") or {}).items():
                if k.startswith("text"): outs.append("".join(v))
        for kind, txt in (("SRC", src), ("OUT", "\n".join(outs))):
            for line in txt.splitlines():
                # report lines with 8-digit numbers 10M..14.9M that are not genomic coordinate rows
                if rx.search(line) and not re.search(r"N[CW]_0\d+\.\d\s", line):
                    print(f"  cell {i} {kind}: {line.strip()[:220]}")
