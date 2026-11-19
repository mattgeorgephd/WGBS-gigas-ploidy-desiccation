#!/usr/bin/env python3
"""Print the bibliographic fields of saved Crossref and Europe PMC records (reader L3).
Usage: python3 -I show_records.py <downloads_dir> <label> [...]
Reads untrusted JSON as data only.
"""
import json
import os
import sys


def main():
    out = sys.argv[1]
    for label in sys.argv[2:]:
        d = os.path.join(out, "records_" + label)
        print("=" * 100)
        print(label)
        try:
            with open(os.path.join(d, "crossref.json"), "rb") as f:
                j = json.loads(f.read().decode("utf-8", "replace"))
            m = j.get("message", {})
            authors = []
            for a in m.get("author", []):
                authors.append((a.get("family") or a.get("name") or "?") + ", " + (a.get("given") or ""))
            print("  CR DOI      :", m.get("DOI"))
            print("  CR type     :", m.get("type"), "| subtype:", m.get("subtype"))
            print("  CR title    :", " | ".join(m.get("title", [])))
            print("  CR authors  :", "; ".join(authors))
            print("  CR container:", " | ".join(m.get("container-title", [])), "| short:", " | ".join(m.get("short-container-title", [])))
            print("  CR volume   :", m.get("volume"), "issue:", m.get("issue"), "page:", m.get("page"), "article-number:", m.get("article-number"))
            for k in ("published-print", "published-online", "published", "issued", "posted"):
                if k in m:
                    print("  CR %-16s:" % k, m[k].get("date-parts"))
            print("  CR publisher:", m.get("publisher"))
            if "relation" in m:
                print("  CR relation :", json.dumps(m["relation"])[:600])
            if "update-to" in m:
                print("  CR update-to:", json.dumps(m["update-to"])[:600])
            if "abstract" in m:
                print("  CR abstract :", m["abstract"][:6000])
            lic = [x.get("URL") for x in m.get("license", [])]
            print("  CR license  :", lic)
            links = [(x.get("URL"), x.get("content-type")) for x in m.get("link", [])]
            print("  CR links    :", links)
        except Exception as e:  # noqa
            print("  crossref: could not parse:", repr(e))
            try:
                with open(os.path.join(d, "crossref.json"), "rb") as f:
                    print("  raw:", f.read()[:300])
            except Exception:
                pass
        try:
            with open(os.path.join(d, "epmc.json"), "rb") as f:
                j = json.loads(f.read().decode("utf-8", "replace"))
            res = j.get("resultList", {}).get("result", [])
            print("  EPMC hits   :", j.get("hitCount"))
            for r in res:
                print("  EPMC id     :", r.get("source"), r.get("id"), "PMID", r.get("pmid"), "PMCID", r.get("pmcid"))
                print("  EPMC title  :", r.get("title"))
                print("  EPMC authors:", r.get("authorString"))
                ji = r.get("journalInfo", {})
                print("  EPMC journal:", (ji.get("journal") or {}).get("title"), "| iso:", (ji.get("journal") or {}).get("isoabbreviation"), "| vol", ji.get("volume"), "issue", ji.get("issue"), "| year", ji.get("yearOfPublication"), "| pages", r.get("pageInfo"))
                print("  EPMC OA     :", r.get("isOpenAccess"), "| inEPMC:", r.get("inEPMC"), "| inPMC:", r.get("inPMC"), "| hasPDF:", r.get("hasPDF"), "| license:", r.get("license"))
                print("  EPMC pubtype:", (r.get("pubTypeList") or {}).get("pubType"))
                ft = (r.get("fullTextUrlList") or {}).get("fullTextUrl", [])
                for u in ft:
                    print("  EPMC url    :", u.get("availability"), u.get("documentStyle"), u.get("site"), u.get("url"))
                print("  EPMC abstract:", r.get("abstractText"))
        except Exception as e:  # noqa
            print("  epmc: could not parse:", repr(e))


if __name__ == "__main__":
    main()
