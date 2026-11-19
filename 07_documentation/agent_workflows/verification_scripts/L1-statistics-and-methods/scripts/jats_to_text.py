# Convert a JATS XML file (Europe PMC fullTextXML) to plain text. No entity expansion beyond the parser defaults; no network.
import sys, re
import xml.etree.ElementTree as ET

src, dst = sys.argv[1], sys.argv[2]
with open(src, "rb") as fh:
    data = fh.read()
# drop a DOCTYPE if present so that the parser does not look for a DTD
data = re.sub(rb"<!DOCTYPE[^>]*>", b"", data, count=1)
root = ET.fromstring(data)

BLOCK = {"p", "title", "sec", "abstract", "caption", "table-wrap", "fig", "ref", "tr", "article-title", "label", "list-item", "disp-formula", "fn", "th", "td", "table-wrap-foot", "kwd", "contrib", "aff", "journal-title", "pub-date", "volume", "issue", "elocation-id", "fpage", "lpage", "article-id", "ack", "notes", "supplementary-material", "def-item", "boxed-text", "statement"}

out = []
def walk(el, depth=0):
    tag = el.tag.split('}')[-1] if isinstance(el.tag, str) else ""
    if tag in ("title",):
        out.append("\n\n## ")
    elif tag in ("td", "th"):
        out.append(" | ")
    elif tag == "tr":
        out.append("\n")
    elif tag in BLOCK:
        out.append("\n")
    if tag == "article-id":
        out.append("[article-id %s] " % el.get("pub-id-type"))
    if tag == "xref" and el.get("ref-type") == "bibr":
        out.append("[")
    if tag == "sup":
        out.append("^")
    if tag == "sub":
        out.append("_")
    if el.text:
        out.append(el.text)
    for ch in el:
        walk(ch, depth + 1)
        if ch.tail:
            out.append(ch.tail)
    if tag == "xref" and el.get("ref-type") == "bibr":
        out.append("]")
    if tag in ("title",):
        out.append("\n")
    elif tag in BLOCK and tag not in ("td", "th"):
        out.append("\n")
walk(root)
txt = "".join(out)
txt = re.sub(r"[ \t]+", " ", txt)
txt = re.sub(r"\n\s*\n\s*\n+", "\n\n", txt)
with open(dst, "w", encoding="utf-8") as fh:
    fh.write(txt)
print(dst, len(txt), "chars")
