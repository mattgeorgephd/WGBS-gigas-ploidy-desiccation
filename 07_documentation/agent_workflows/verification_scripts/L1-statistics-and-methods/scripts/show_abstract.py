import json, sys, re
for fn in sys.argv[1:]:
    j = json.load(open(fn, encoding="utf-8"))
    for r in j.get("resultList", {}).get("result", []):
        print("==", r.get("authorString"), r.get("pubYear"), r.get("title"))
        print(re.sub(r"<[^>]+>", "", r.get("abstractText") or "(no abstract)"))
        print()
