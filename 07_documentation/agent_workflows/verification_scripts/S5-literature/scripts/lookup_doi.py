# Look up published papers by DOI at Crossref and Europe PMC.
# Sends only the DOI. Saves the raw JSON (untrusted data, parsed as JSON only)
# under the output directory given as argv[1]; prints the metadata fields.
# Usage: python3 -I lookup_doi.py <outdir> <label>=<doi> [<label>=<doi> ...]
import json
import os
import subprocess
import sys
import time
import urllib.parse

UA = "literature-metadata-check/1.0"


def fetch(url, dest, tries=3, pause=6):
    code = "000"
    for i in range(tries):
        r = subprocess.run(
            ["curl", "-sS", "-m", "60", "-A", UA, "-o", dest, "-w", "%{http_code}", url],
            capture_output=True, text=True)
        code = r.stdout.strip()
        if code == "200":
            return code
        if code in ("403", "407"):
            return code  # policy or refusal: do not retry
        time.sleep(pause * (i + 1))
    return code


def show_crossref(path):
    try:
        d = json.load(open(path, encoding="utf-8"))["message"]
    except Exception as e:
        print("  crossref: not JSON (%s)" % e)
        return
    auth = "; ".join(
        ("%s %s" % (a.get("family", a.get("name", "")), a.get("given", ""))).strip()
        for a in d.get("author", []))
    print("  CR title   :", " | ".join(d.get("title", [])))
    print("  CR authors :", auth)
    print("  CR journal :", " | ".join(d.get("container-title", [])),
          "| short:", " | ".join(d.get("short-container-title", [])))
    print("  CR vol/iss/page/artno:", d.get("volume"), "/", d.get("issue"), "/",
          d.get("page"), "/", d.get("article-number"))
    for k in ("published-print", "published-online", "published", "issued"):
        if d.get(k):
            print("  CR %-16s:" % k, d[k].get("date-parts"))
    print("  CR type/publisher:", d.get("type"), "/", d.get("publisher"))
    print("  CR DOI     :", d.get("DOI"))
    if d.get("relation"):
        print("  CR relation:", json.dumps(d.get("relation"))[:400])
    if d.get("update-to") or d.get("updated-by"):
        print("  CR updates :", json.dumps(d.get("update-to")), json.dumps(d.get("updated-by")))
    if d.get("abstract"):
        print("  CR abstract:", d["abstract"][:3000])


def show_epmc(path):
    try:
        d = json.load(open(path, encoding="utf-8"))
    except Exception as e:
        print("  epmc: not JSON (%s)" % e)
        return
    print("  EPMC hits:", d.get("hitCount"))
    for r in d.get("resultList", {}).get("result", []):
        ji = r.get("journalInfo") or {}
        j = ji.get("journal") or {}
        print("  EPMC id/source/pmid/pmcid:", r.get("id"), r.get("source"), r.get("pmid"), r.get("pmcid"))
        print("  EPMC title  :", r.get("title"))
        print("  EPMC authors:", r.get("authorString"))
        print("  EPMC journal:", j.get("isoabbreviation"), "| vol", ji.get("volume"), "| issue", ji.get("issue"),
              "| pages", r.get("pageInfo"), "| year", r.get("pubYear"), "| date", ji.get("dateOfPublication"))
        print("  EPMC OA/inEPMC/license:", r.get("isOpenAccess"), r.get("inEPMC"), r.get("license"))
        print("  EPMC abstract:", r.get("abstractText"))


def main():
    outdir = sys.argv[1]
    os.makedirs(os.path.join(outdir, "crossref"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "epmc"), exist_ok=True)
    for item in sys.argv[2:]:
        label, doi = item.split("=", 1)
        print("=" * 100)
        print(label, doi)
        cr = os.path.join(outdir, "crossref", label + ".json")
        code = fetch("https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="/"), cr)
        print("  crossref http", code)
        if code == "200":
            show_crossref(cr)
        time.sleep(2)
        ep = os.path.join(outdir, "epmc", label + ".json")
        q = urllib.parse.quote('DOI:"%s"' % doi, safe="")
        code = fetch("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=" + q +
                     "&resultType=core&format=json", ep)
        print("  epmc http", code)
        if code == "200":
            show_epmc(ep)
        time.sleep(2)


if __name__ == "__main__":
    main()
