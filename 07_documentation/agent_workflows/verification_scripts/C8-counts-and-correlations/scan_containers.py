"""Search every file under given roots for the manuscript's CpG count and correlation range.
Zip containers (xlsx, pptx, docx) are opened and every member searched; nested zips are recursed.
PDFs are converted with pdftotext. Everything else is searched as raw bytes (latin-1)."""
import sys, os, re, zipfile, io, subprocess, collections

EXACT = re.compile(r"12[\s,.'  ]?794[\s,.'  ]?074|1\.2794074E\+?0?7|1\.2794074e\+?0?7")
# nearby: integers 12,600,000..12,999,999 written with or without separators; and 12.7x / 12.8x million
NEAR = re.compile(r"(?<![\d.,])12[\s,]?[6-9]\d{2}[\s,]?\d{3}(?![\d])")
MILLION = re.compile(r"12\.[6-9]\d?\s?(?:M\b|million|Mio)", re.I)
PCT = re.compile(r"96\.[56]\d*\s?%?")
COR = re.compile(r"0\.87\b.{0,40}0\.93\b|0\.93\b.{0,40}0\.87\b", re.S)
COR1 = re.compile(r"(?<![\d.])0?\.(87|93)(?![\d])")

stats = collections.Counter()
hits = []

def search_text(label, text):
    stats["texts"] += 1
    for name, rx in (("EXACT", EXACT), ("NEAR", NEAR), ("MILLION", MILLION), ("COR_RANGE", COR)):
        for m in rx.finditer(text):
            s = max(0, m.start() - 60); e = min(len(text), m.end() + 60)
            hits.append((name, label, text[s:e].replace("\n", " ")))
            if name != "EXACT" and sum(1 for h in hits if h[0] == name and h[1] == label) > 8:
                break

def handle_bytes(label, data, depth=0):
    if data[:4] == b"PK\x03\x04" and depth < 4:
        try:
            z = zipfile.ZipFile(io.BytesIO(data))
        except Exception as ex:
            stats["bad_zip"] += 1; return
        stats["zip_containers"] += 1
        for info in z.infolist():
            try:
                d = z.read(info)
            except Exception:
                stats["bad_member"] += 1; continue
            stats["zip_members"] += 1
            handle_bytes(label + "!" + info.filename, d, depth + 1)
        return
    if data[:5] == b"%PDF-":
        stats["pdf"] += 1
        try:
            out = subprocess.run(["pdftotext", "-layout", "-", "-"], input=data, capture_output=True, timeout=120).stdout
            search_text(label + " [pdftotext]", out.decode("utf-8", "replace"))
        except Exception as ex:
            stats["pdf_fail"] += 1
        return
    # raw search; also try utf-16 for old .xls (BIFF stores some strings as utf-16le)
    txt = data.decode("latin-1")
    search_text(label, txt)
    if label.lower().endswith((".xls", ".ppt", ".doc")):
        search_text(label + " [utf16]", data.decode("utf-16-le", "replace"))

def walk(root):
    for dp, dn, fn in os.walk(root):
        if ".git" in dp.split(os.sep): continue
        for f in fn:
            p = os.path.join(dp, f)
            ext = os.path.splitext(f)[1].lower()
            stats["files"] += 1; stats["ext" + ext] += 1
            if ext in (".jpeg", ".jpg", ".png", ".tiff", ".tif"):
                stats["images_not_searchable"] += 1
                continue
            try:
                with open(p, "rb") as fh: data = fh.read()
            except Exception as ex:
                stats["unreadable"] += 1; continue
            handle_bytes(p, data)

for r in sys.argv[1:]:
    walk(r)
print("STATS", dict(stats))
by = collections.defaultdict(list)
for h in hits: by[h[0]].append(h)
for k in ("EXACT", "NEAR", "MILLION", "COR_RANGE"):
    print(f"\n==== {k}: {len(by[k])} hits")
    seen = collections.Counter()
    for _, label, ctx in by[k]:
        seen[label] += 1
        if seen[label] <= 3:
            print("  ", label[-110:], "::", ctx[:200])
