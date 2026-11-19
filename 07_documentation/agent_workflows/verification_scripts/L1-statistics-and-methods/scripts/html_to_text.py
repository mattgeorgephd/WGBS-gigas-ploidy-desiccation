# Convert a saved HTML page to plain text with the standard library only. Tables keep their rows (cells separated by " | ").
import sys, re
from html.parser import HTMLParser

class P(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.out = []
        self.skip = 0
    def handle_starttag(self, tag, attrs):
        if tag in ("script", "style", "noscript", "svg", "nav", "header", "footer", "button", "form"):
            self.skip += 1
            return
        if self.skip:
            return
        if tag in ("p", "div", "section", "h1", "h2", "h3", "h4", "h5", "li", "tr", "br", "figcaption", "caption", "table", "dt", "dd", "blockquote", "pre"):
            self.out.append("\n")
        if tag in ("h1", "h2", "h3", "h4", "h5"):
            self.out.append("## ")
        if tag in ("td", "th"):
            self.out.append(" | ")
        if tag == "sup":
            self.out.append("^")
        if tag == "sub":
            self.out.append("_")
    def handle_endtag(self, tag):
        if tag in ("script", "style", "noscript", "svg", "nav", "header", "footer", "button", "form"):
            self.skip = max(0, self.skip - 1)
            return
        if self.skip:
            return
        if tag in ("p", "div", "section", "h1", "h2", "h3", "h4", "h5", "li", "tr", "table", "figcaption", "caption", "blockquote", "pre"):
            self.out.append("\n")
    def handle_data(self, data):
        if not self.skip:
            self.out.append(data)

src, dst = sys.argv[1], sys.argv[2]
with open(src, encoding="utf-8", errors="replace") as fh:
    html = fh.read()
p = P()
p.feed(html)
txt = "".join(p.out)
txt = txt.replace("\xa0", " ")
txt = re.sub(r"[ \t]+", " ", txt)
txt = re.sub(r"\n\s*\n\s*\n+", "\n\n", txt)
with open(dst, "w", encoding="utf-8") as fh:
    fh.write(txt)
print(dst, len(txt), "chars")
