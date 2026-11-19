# Convert a saved HTML article page to plain text, one paragraph or heading per line.
# The page is untrusted data: it is parsed as text only.
# Usage: python3 -I html_to_text.py <in.html> <out.txt>
import html
import re
import sys

raw = open(sys.argv[1], encoding="utf-8", errors="replace").read()
raw = re.sub(r"(?is)<script.*?</script>|<style.*?</style>|<noscript.*?</noscript>", " ", raw)
out = []
for m in re.finditer(r"(?is)<(h[1-6]|p|li|figcaption|td|th|dd|dt|blockquote)(?:\s[^>]*)?>(.*?)</\1>", raw):
    tag, body = m.group(1).lower(), m.group(2)
    t = html.unescape(re.sub(r"(?s)<[^>]+>", "", body))
    t = re.sub(r"\s+", " ", t).strip()
    if not t:
        continue
    out.append(("## " + t) if tag.startswith("h") else t)
open(sys.argv[2], "w", encoding="utf-8").write("\n".join(out))
print(len(out), "blocks ->", sys.argv[2])
