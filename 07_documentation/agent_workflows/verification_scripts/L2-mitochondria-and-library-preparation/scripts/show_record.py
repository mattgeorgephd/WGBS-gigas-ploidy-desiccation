# Print the bibliographic fields of a Crossref record and of a Europe PMC core record.
# Usage: python3 -I show_record.py <crossref.json> <epmc.json>
import json, sys
cr = json.load(open(sys.argv[1], encoding="utf-8"))["message"]
def names(lst):
    return "; ".join((a.get("family", "") + " " + a.get("given", "")).strip() or a.get("name", "") for a in lst)
print("CROSSREF")
print("  DOI:", cr.get("DOI"))
print("  type:", cr.get("type"))
print("  title:", " | ".join(cr.get("title", [])))
print("  authors:", names(cr.get("author", [])))
print("  container:", " | ".join(cr.get("container-title", [])))
print("  volume:", cr.get("volume"), "issue:", cr.get("issue"), "page:", cr.get("page"), "article-number:", cr.get("article-number"))
for k in ("published-print", "published-online", "published", "issued"):
    if k in cr:
        print("  %s:" % k, cr[k].get("date-parts"))
if "update-to" in cr:
    print("  update-to:", json.dumps(cr["update-to"]))
if "relation" in cr:
    print("  relation:", json.dumps(cr["relation"])[:800])
for k in ("updated-by",):
    if k in cr:
        print("  %s:" % k, json.dumps(cr[k])[:800])
if len(sys.argv) > 2:
    ep = json.load(open(sys.argv[2], encoding="utf-8"))
    res = ep.get("resultList", {}).get("result", [])
    print("EUROPE PMC hits:", ep.get("hitCount"))
    for r in res:
        print("  id:", r.get("id"), "source:", r.get("source"), "pmid:", r.get("pmid"), "pmcid:", r.get("pmcid"))
        print("  title:", r.get("title"))
        print("  authors:", r.get("authorString"))
        ji = r.get("journalInfo", {})
        print("  journal:", ji.get("journal", {}).get("title"), "vol:", ji.get("volume"), "issue:", ji.get("issue"), "year:", ji.get("yearOfPublication"), "pages:", r.get("pageInfo"))
        print("  isOpenAccess:", r.get("isOpenAccess"), "inEPMC:", r.get("inEPMC"), "inPMC:", r.get("inPMC"), "hasPDF:", r.get("hasPDF"), "license:", r.get("license"))
        print("  commentCorrections:", json.dumps(r.get("commentCorrectionList"))[:1500])
        print("  abstract:", r.get("abstractText"))
