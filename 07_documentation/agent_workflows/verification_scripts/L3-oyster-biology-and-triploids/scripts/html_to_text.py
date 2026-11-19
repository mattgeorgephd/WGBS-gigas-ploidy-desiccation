#!/usr/bin/env python3
"""Convert a saved HTML page to plain text (reader L3). The page is untrusted data:
it is parsed with the standard-library HTML parser only; scripts and styles are dropped.
Usage: python3 -I html_to_text.py <in.html> <out.txt>
"""
import re
import sys
from html.parser import HTMLParser

SKIP = {"script", "style", "noscript", "svg", "template"}
BLOCK = {"p", "div", "section", "h1", "h2", "h3", "h4", "h5", "li", "tr", "br", "table", "figcaption",
         "article", "header", "footer", "dt", "dd", "ul", "ol", "blockquote", "title"}


class P(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.out = []
        self.skip = 0
        self.meta = []

    def handle_starttag(self, tag, attrs):
        if tag in SKIP:
            self.skip += 1
        if tag in BLOCK:
            self.out.append("\n")
        if tag in ("h1", "h2", "h3", "h4"):
            self.out.append("## ")
        if tag in ("td", "th"):
            self.out.append(" | ")
        if tag == "meta":
            d = dict(attrs)
            n = d.get("name") or d.get("property")
            if n and d.get("content") and re.match(r"(citation_|dc\.|DC\.|description|prism\.|og:title|og:description)", n):
                self.meta.append("%s = %s" % (n, d.get("content")))

    def handle_endtag(self, tag):
        if tag in SKIP and self.skip:
            self.skip -= 1
        if tag in BLOCK:
            self.out.append("\n")

    def handle_data(self, data):
        if not self.skip:
            self.out.append(data)


def main():
    with open(sys.argv[1], "rb") as f:
        raw = f.read().decode("utf-8", "replace")
    p = P()
    p.feed(raw)
    text = "".join(p.out)
    text = re.sub(r"[ \t\r\f\v]+", " ", text)
    text = re.sub(r"\n[ \t]+", "\n", text)
    text = re.sub(r"\n{3,}", "\n\n", text)
    with open(sys.argv[2], "w", encoding="utf-8") as f:
        f.write("META\n" + "\n".join(p.meta) + "\n\nTEXT\n" + text)
    print(sys.argv[2], len(text), "chars;", len(p.meta), "meta")


if __name__ == "__main__":
    main()
