# Convert a JATS XML file to plain text without an XML parser (tags are stripped by regular expressions,
# no entity expansion beyond the standard character references). Usage: python3 -I jats_to_text.py <in.xml> <out.txt>
import re, sys, html
src = open(sys.argv[1], encoding="utf-8", errors="replace").read()
m = re.search(r"<ref-list.*?</ref-list>", src, flags=re.S)
refs = m.group(0) if m else ""
body = src.replace(refs, "") if refs else src
BLOCK = r"p|sec|caption|table-wrap|fig|tr|abstract|ref|list-item|disp-formula|fn|contrib|aff|article-title|kwd|supplementary-material|ack|notes"
def clean(x):
    x = re.sub(r"<xref\b[^>]*ref-type=\"bibr\"[^>]*>(.*?)</xref>", r"[\1]", x, flags=re.S)
    x = re.sub(r"<(?:title|label)(?:\s[^>]*)?>", "\n\n## ", x)
    x = re.sub(r"</(?:title|label)>", "\n", x)
    x = re.sub(r"<(?:%s)(?:\s[^>]*)?>" % BLOCK, "\n", x)
    x = re.sub(r"</(?:%s)>" % BLOCK, "\n", x)
    x = re.sub(r"</(?:td|th)>", " | ", x)
    x = re.sub(r"<[^>]+>", "", x)
    x = html.unescape(x)
    x = re.sub(r"[ \t]+", " ", x)
    x = re.sub(r"\n\s*\n\s*\n+", "\n\n", x)
    return x.strip()
out = clean(body) + "\n\n#### REFERENCES ####\n" + clean(refs)
open(sys.argv[2], "w", encoding="utf-8").write(out)
print(sys.argv[2], len(out), "characters")
