import json, sys
for fn in sys.argv[1:]:
    j = json.load(open(fn, encoding="utf-8"))
    print("==", fn.split("/")[-1], "hits:", j.get("hitCount"))
    for r in j.get("resultList", {}).get("result", []):
        print("  %s | %s | %s %s | %s | doi:%s | pmcid:%s | OA:%s" % (r.get("pubYear"), r.get("authorString", "")[:60], r.get("journalTitle"), r.get("journalVolume"), r.get("title", "")[:150], r.get("doi"), r.get("pmcid"), r.get("isOpenAccess")))
