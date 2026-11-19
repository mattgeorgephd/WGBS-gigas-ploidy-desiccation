# Print open-access locations from an OpenAlex work record. Usage: python3 -I show_openalex.py <openalex.json>
import json, sys
w = json.load(open(sys.argv[1], encoding="utf-8"))
print("title:", w.get("title"))
print("year:", w.get("publication_year"), "| cited_by:", w.get("cited_by_count"))
oa = w.get("open_access", {})
print("open_access:", oa)
for loc in w.get("locations", []):
    src = (loc.get("source") or {}).get("display_name")
    print(" -", src, "| is_oa:", loc.get("is_oa"), "| version:", loc.get("version"), "| landing:", loc.get("landing_page_url"), "| pdf:", loc.get("pdf_url"))
