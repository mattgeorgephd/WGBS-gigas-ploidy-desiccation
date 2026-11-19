import subprocess, sys, re, zipfile, io, collections
EXACT = re.compile(rb"12[\s,.']?794[\s,.']?074|1\.2794074[eE]\+?0?7")
COR = re.compile(rb"0\.87\b.{0,40}0\.93\b|0\.93\b.{0,40}0\.87\b", re.S)
repo = sys.argv[1]; objs = sys.argv[2]
stats = collections.Counter(); hits = []
seen = set()
def scan(label, data, depth=0):
    if data[:4] == b"PK\x03\x04" and depth < 4:
        try: z = zipfile.ZipFile(io.BytesIO(data))
        except Exception: stats["badzip"] += 1; return
        stats["zips"] += 1
        for info in z.infolist():
            try: d = z.read(info)
            except Exception: continue
            scan(label + "!" + info.filename, d, depth + 1)
        return
    if data[:5] == b"%PDF-":
        try:
            data = subprocess.run(["pdftotext", "-layout", "-", "-"], input=data, capture_output=True, timeout=120).stdout
        except Exception: pass
    for name, rx in (("EXACT", EXACT), ("COR", COR)):
        for m in rx.finditer(data):
            hits.append((name, label, data[max(0, m.start()-80):m.end()+80]))
for line in open(objs, encoding="utf-8", errors="replace"):
    sha, path = line.rstrip("\n").split(" ", 1)
    if sha in seen: continue
    seen.add(sha)
    t = subprocess.run(["git", "-C", repo, "cat-file", "-t", sha], capture_output=True).stdout.strip()
    if t != b"blob": continue
    low = path.lower()
    if low.endswith((".jpeg", ".jpg", ".png", ".tiff", ".tif")): stats["image_blobs_skipped"] += 1; continue
    size = int(subprocess.run(["git", "-C", repo, "cat-file", "-s", sha], capture_output=True).stdout.strip())
    if size > 300_000_000: stats["too_big"] += 1; continue
    data = subprocess.run(["git", "-C", repo, "cat-file", "blob", sha], capture_output=True).stdout
    stats["blobs_scanned"] += 1
    scan(f"{sha[:8]}:{path}", data)
print(dict(stats))
audit_prefixes = ("issues.md", "05_manuscript/", "07_documentation/", "03_input_files/manuscript/", "06_diagnostics/", "04_output/", "01_pipeline/", "README.md")
for name in ("EXACT", "COR"):
    hs = [h for h in hits if h[0] == name]
    print(f"== {name}: {len(hs)} hits")
    labs = collections.Counter(h[1] for h in hs)
    for lab, n in labs.most_common(60):
        p = lab.split(":", 1)[1]
        tag = "AUDIT-ERA" if p.startswith(audit_prefixes) else "LEGACY?"
        print(f"   {tag} {n:3d} {lab[:150]}")
    for h in hs:
        p = h[1].split(":", 1)[1]
        if not p.startswith(audit_prefixes):
            print("   CONTEXT", h[1][:120], "::", h[2][:240])
