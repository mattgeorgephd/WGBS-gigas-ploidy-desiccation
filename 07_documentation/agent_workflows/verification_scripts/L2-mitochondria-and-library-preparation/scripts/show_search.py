# List hits of a Europe PMC search result. Usage: python3 -I show_search.py <search.json>
import json, sys
r = json.load(open(sys.argv[1], encoding="utf-8"))
print("hitCount:", r.get("hitCount"))
for x in r.get("resultList", {}).get("result", []):
    print("-", x.get("pubYear"), "|", x.get("authorString", "")[:60], "|", x.get("title", "")[:170], "|", x.get("journalTitle"), x.get("journalVolume"), x.get("pageInfo"), "| doi:", x.get("doi"), "| pmcid:", x.get("pmcid"), "| OA:", x.get("isOpenAccess"), "| cited:", x.get("citedByCount"))
