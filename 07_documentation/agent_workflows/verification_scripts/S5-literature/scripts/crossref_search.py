# Bibliographic search at Crossref with generic scientific terms; lists the hits.
# Sends only the query string. Saves raw JSON (untrusted, parsed as JSON only).
# Usage: python3 -I crossref_search.py <outdir> <label> <query> [rows] [extra filter]
import json
import os
import subprocess
import sys
import time
import urllib.parse

UA = "literature-metadata-check/1.0"
outdir, label, query = sys.argv[1], sys.argv[2], sys.argv[3]
rows = sys.argv[4] if len(sys.argv) > 4 else "20"
flt = sys.argv[5] if len(sys.argv) > 5 else ""
os.makedirs(os.path.join(outdir, "crossref_search"), exist_ok=True)
dest = os.path.join(outdir, "crossref_search", label + ".json")
url = ("https://api.crossref.org/works?query.bibliographic=" + urllib.parse.quote(query, safe="") +
       "&rows=" + rows + "&select=DOI,title,author,container-title,volume,page,article-number,issued,type")
if flt:
    url += "&filter=" + urllib.parse.quote(flt, safe=":,-")
code = "000"
for i in range(3):
    r = subprocess.run(["curl", "-sS", "-m", "60", "-A", UA, "-o", dest, "-w", "%{http_code}", url],
                       capture_output=True, text=True)
    code = r.stdout.strip()
    if code == "200" or code in ("403", "407"):
        break
    time.sleep(8 * (i + 1))
print("query:", query, "| http", code)
if code == "200":
    d = json.load(open(dest, encoding="utf-8"))["message"]
    for it in d.get("items", []):
        au = it.get("author", [])
        a = (au[0].get("family", au[0].get("name", "")) if au else "") + (" et al." if len(au) > 1 else "")
        yr = (it.get("issued", {}).get("date-parts") or [[None]])[0][0]
        print("- %s %s | %s | %s %s:%s | doi:%s | %s" % (
            a, yr, " ".join(it.get("title", []))[:170], " ".join(it.get("container-title", []))[:50],
            it.get("volume", ""), it.get("page", it.get("article-number", "")), it.get("DOI", ""), it.get("type", "")))
