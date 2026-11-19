import json, sys, os
d = sys.argv[1]
for fn in sorted(os.listdir(d)):
    if not fn.endswith(".json"):
        continue
    with open(os.path.join(d, fn), encoding="utf-8") as fh:
        j = json.load(fh)
    res = j.get("resultList", {}).get("result", [])
    print("==", fn, "hits:", j.get("hitCount"))
    for r in res:
        print("  id:", r.get("id"), "source:", r.get("source"), "pmid:", r.get("pmid"), "pmcid:", r.get("pmcid"), "doi:", r.get("doi"))
        print("  title:", r.get("title"))
        print("  authors:", r.get("authorString"))
        ji = r.get("journalInfo", {})
        print("  journal:", ji.get("journal", {}).get("title"), "vol:", ji.get("volume"), "issue:", ji.get("issue"), "year:", ji.get("yearOfPublication"), "pageInfo:", r.get("pageInfo"))
        print("  isOpenAccess:", r.get("isOpenAccess"), "inEPMC:", r.get("inEPMC"), "inPMC:", r.get("inPMC"), "hasPDF:", r.get("hasPDF"), "license:", r.get("license"))
        ab = r.get("abstractText") or ""
        print("  abstract chars:", len(ab))
        ftl = r.get("fullTextUrlList", {}).get("fullTextUrl", [])
        for u in ftl:
            print("    url:", u.get("availability"), u.get("documentStyle"), u.get("site"), u.get("url"))
