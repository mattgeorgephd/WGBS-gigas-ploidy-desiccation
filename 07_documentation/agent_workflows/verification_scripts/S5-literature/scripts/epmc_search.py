# Search Europe PMC by a query string of generic scientific terms and list the hits.
# Sends only the query. Saves raw JSON (untrusted, parsed as JSON only).
# Usage: python3 -I epmc_search.py <outdir> <label> <query> [pageSize]
import json
import os
import subprocess
import sys
import time
import urllib.parse

UA = "literature-metadata-check/1.0"
outdir, label, query = sys.argv[1], sys.argv[2], sys.argv[3]
n = sys.argv[4] if len(sys.argv) > 4 else "25"
os.makedirs(os.path.join(outdir, "epmc_search"), exist_ok=True)
dest = os.path.join(outdir, "epmc_search", label + ".json")
url = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=" + urllib.parse.quote(query, safe="") +
       "&resultType=lite&format=json&pageSize=" + n)
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
    d = json.load(open(dest, encoding="utf-8"))
    print("hits:", d.get("hitCount"))
    for r in d.get("resultList", {}).get("result", []):
        print("- %s | %s | %s %s;%s:%s | doi:%s | %s%s | OA:%s" % (
            r.get("authorString", "")[:70], r.get("title", "")[:170], r.get("journalTitle", r.get("bookOrReportDetails", "")),
            r.get("pubYear", ""), r.get("journalVolume", ""), r.get("pageInfo", ""), r.get("doi", ""),
            r.get("source", ""), r.get("id", ""), r.get("isOpenAccess", "")))
