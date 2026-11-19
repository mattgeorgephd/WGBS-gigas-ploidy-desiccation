#!/usr/bin/env python3
"""Literature searches with generic terms (reader L3): Crossref and Europe PMC.
Usage:
  python3 -I search.py crossref <downloads_dir> <label> <rows> <query words...>
  python3 -I search.py epmc     <downloads_dir> <label> <rows> <query string>
Only the query string is sent. Generic client name. One request per call.
"""
import json
import os
import sys
import urllib.parse
import urllib.request

UA = "literature-check/1.0"


def get(url, timeout=90):
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.status, r.read()
    except urllib.error.HTTPError as e:
        return e.code, e.read()
    except Exception as e:  # noqa
        return -1, repr(e).encode()


def main():
    mode, dl, label, rows = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    query = " ".join(sys.argv[5:])
    d = os.path.join(dl, "search_" + mode + "_" + label)
    os.makedirs(d, exist_ok=True)
    if mode == "crossref":
        url = "https://api.crossref.org/works?rows=%d&query.bibliographic=%s&select=DOI,title,author,container-title,issued,type,volume,page" % (
            rows, urllib.parse.quote(query))
    else:
        url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=%s&format=json&resultType=lite&pageSize=%d" % (
            urllib.parse.quote(query), rows)
    st, body = get(url)
    with open(os.path.join(d, "result.json"), "wb") as f:
        f.write(body)
    print("query:", query, "| http", st)
    try:
        j = json.loads(body.decode("utf-8", "replace"))
    except Exception as e:  # noqa
        print("  parse error", repr(e), body[:200])
        return
    if mode == "crossref":
        for it in j.get("message", {}).get("items", []):
            au = it.get("author", [])
            first = (au[0].get("family") or au[0].get("name") or "?") if au else "?"
            yr = (it.get("issued", {}).get("date-parts") or [[None]])[0][0]
            print("- %s et al. %s | %s | %s %s:%s | doi:%s | %s" % (
                first, yr, " ".join(it.get("title", [""]))[:170], " ".join(it.get("container-title", [""])),
                it.get("volume", ""), it.get("page", ""), it.get("DOI"), it.get("type")))
    else:
        print("hits:", j.get("hitCount"))
        for r in j.get("resultList", {}).get("result", []):
            print("- %s | %s | %s %s;%s:%s | doi:%s | %s%s | OA:%s | PMCID:%s" % (
                (r.get("authorString") or "")[:70], (r.get("title") or "")[:180], r.get("journalTitle") or r.get("bookOrReportDetails"),
                r.get("pubYear"), r.get("journalVolume"), r.get("pageInfo"), r.get("doi"), r.get("source"), r.get("id"),
                r.get("isOpenAccess"), r.get("pmcid")))


if __name__ == "__main__":
    main()
