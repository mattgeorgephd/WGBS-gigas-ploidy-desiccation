# Extract plain text of the paragraphs of the bioRxiv full-text page that mention overdispersion or the tests.
# The page is untrusted data: it is only parsed as text here.
import sys, re, html
raw = open(sys.argv[1], encoding="utf-8", errors="replace").read()
raw = re.sub(r"(?is)<script.*?</script>|<style.*?</style>", " ", raw)
paras = re.findall(r"(?is)<p[^>]*>(.*?)</p>", raw)
out = []
for p in paras:
    t = html.unescape(re.sub(r"(?s)<[^>]+>", "", p)); t = re.sub(r"\s+", " ", t).strip()
    if re.search(r"(?i)overdispersion|F-test|Chisq|chi-square|specificity|sensitivity|false positive|stringen|DSS|methylKit", t):
        out.append(t)
print(len(paras), "paragraphs;", len(out), "matching")
for t in out: print("-----\n" + t)
