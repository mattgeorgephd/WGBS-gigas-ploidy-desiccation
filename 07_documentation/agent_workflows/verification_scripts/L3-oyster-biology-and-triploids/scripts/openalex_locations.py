#!/usr/bin/env python3
"""List open-access locations of a work from OpenAlex by DOI (reader L3).
Usage: python3 -I openalex_locations.py <downloads_dir> <label>=<doi> [...]
Only the DOI is sent; generic client name; one request at a time.
"""
import json
import os
import sys
import time
import urllib.parse
import urllib.request

UA = "literature-check/1.0"


def get(url, timeout=60):
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.status, r.read()
    except urllib.error.HTTPError as e:
        return e.code, e.read()
    except Exception as e:  # noqa
        return -1, repr(e).encode()


def main():
    dl = sys.argv[1]
    for item in sys.argv[2:]:
        label, doi = item.split("=", 1)
        d = os.path.join(dl, "openalex_" + label)
        os.makedirs(d, exist_ok=True)
        url = "https://api.openalex.org/works/doi:" + urllib.parse.quote(doi, safe="/")
        st, body = get(url)
        with open(os.path.join(d, "work.json"), "wb") as f:
            f.write(body)
        print("=" * 80)
        print(label, doi, "http", st)
        try:
            j = json.loads(body.decode("utf-8", "replace"))
        except Exception as e:  # noqa
            print("  parse error", repr(e), body[:200])
            continue
        print("  title:", j.get("title"))
        print("  year:", j.get("publication_year"), "| type:", j.get("type"), "| cited_by:", j.get("cited_by_count"))
        oa = j.get("open_access", {})
        print("  open_access:", oa.get("is_oa"), oa.get("oa_status"), oa.get("oa_url"))
        for loc in j.get("locations", []) or []:
            src = (loc.get("source") or {})
            print("  location:", src.get("display_name"), "|", src.get("type"), "| oa:", loc.get("is_oa"), "| version:", loc.get("version"),
                  "| landing:", loc.get("landing_page_url"), "| pdf:", loc.get("pdf_url"))
        time.sleep(1.0)


if __name__ == "__main__":
    main()
