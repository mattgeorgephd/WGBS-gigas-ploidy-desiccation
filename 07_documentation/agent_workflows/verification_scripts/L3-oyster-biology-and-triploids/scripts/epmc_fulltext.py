#!/usr/bin/env python3
"""Fetch open-access full text XML from Europe PMC and write plain text (reader L3).
Usage: python3 -I epmc_fulltext.py <downloads_dir> <text_dir> <label>=<PMCID> [...]
One request at a time; generic client name. The XML is parsed as data only
(no DTD loading, no entity expansion beyond the built-in parser defaults).
"""
import os
import re
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET

UA = "literature-check/1.0"


def get(url, timeout=120):
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.status, r.read()
    except urllib.error.HTTPError as e:
        return e.code, e.read()
    except Exception as e:  # noqa
        return -1, repr(e).encode()


BLOCK = {"p", "title", "sec", "caption", "table-wrap", "fig", "tr", "ref", "abstract", "article-title",
         "label", "list-item", "td", "th", "table", "thead", "tbody", "disp-formula", "supplementary-material",
         "kwd", "kwd-group", "contrib", "aff", "fn", "ack", "body", "back", "front", "journal-title", "pub-date",
         "volume", "issue", "elocation-id", "fpage", "lpage", "article-id", "subject", "def-item", "term", "def"}


def localname(tag):
    return tag.split("}")[-1] if isinstance(tag, str) else ""


def walk(el, out):
    name = localname(el.tag)
    if name in ("tr",):
        cells = []
        for c in el:
            cells.append(re.sub(r"\s+", " ", "".join(c.itertext())).strip())
        out.append("\n| " + " | ".join(cells) + " |")
        if el.tail:
            out.append(el.tail)
        return
    if name == "title":
        out.append("\n\n## ")
    elif name in BLOCK:
        out.append("\n")
    if name == "article-id":
        out.append("[article-id %s] " % el.attrib.get("pub-id-type", ""))
    if el.text:
        out.append(el.text)
    for c in el:
        walk(c, out)
    if name in BLOCK:
        out.append("\n")
    if el.tail:
        out.append(el.tail)


def main():
    dl, td = sys.argv[1], sys.argv[2]
    for item in sys.argv[3:]:
        label, pmcid = item.split("=", 1)
        d = os.path.join(dl, "epmc_fulltext_" + label)
        os.makedirs(d, exist_ok=True)
        url = "https://www.ebi.ac.uk/europepmc/webservices/rest/%s/fullTextXML" % pmcid
        st, body = get(url)
        xmlpath = os.path.join(d, pmcid + ".xml")
        with open(xmlpath, "wb") as f:
            f.write(body)
        print(label, pmcid, "status", st, "bytes", len(body))
        if st != 200:
            continue
        try:
            root = ET.fromstring(body)
        except Exception as e:  # noqa
            print("  parse error", repr(e))
            continue
        out = []
        walk(root, out)
        text = "".join(out)
        text = re.sub(r"[ \t]+", " ", text)
        text = re.sub(r"\n[ \t]+", "\n", text)
        text = re.sub(r"\n{3,}", "\n\n", text)
        tp = os.path.join(td, "%s_%s.txt" % (label, pmcid))
        with open(tp, "w", encoding="utf-8") as f:
            f.write("SOURCE: %s\n\n" % url)
            f.write(text)
        print("  text ->", tp, len(text))
        time.sleep(1.5)


if __name__ == "__main__":
    main()
