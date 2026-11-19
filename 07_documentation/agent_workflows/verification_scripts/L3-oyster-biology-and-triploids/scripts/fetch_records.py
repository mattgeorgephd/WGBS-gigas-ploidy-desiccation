#!/usr/bin/env python3
"""Fetch Crossref and Europe PMC records by DOI (reader L3). One request at a time.
Usage: python3 -I fetch_records.py <downloads_dir> <label>=<doi> [...]
Only DOIs are sent. Generic client name.
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
    out = sys.argv[1]
    for item in sys.argv[2:]:
        label, doi = item.split("=", 1)
        d = os.path.join(out, "records_" + label)
        os.makedirs(d, exist_ok=True)
        # Crossref
        url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="/")
        st, body = get(url)
        with open(os.path.join(d, "crossref.json"), "wb") as f:
            f.write(body)
        print(label, "crossref", st, len(body))
        time.sleep(1.0)
        # Europe PMC
        q = urllib.parse.quote('DOI:"%s"' % doi)
        url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=%s&format=json&resultType=core" % q
        st, body = get(url)
        with open(os.path.join(d, "epmc.json"), "wb") as f:
            f.write(body)
        print(label, "epmc", st, len(body))
        time.sleep(1.0)


if __name__ == "__main__":
    main()
