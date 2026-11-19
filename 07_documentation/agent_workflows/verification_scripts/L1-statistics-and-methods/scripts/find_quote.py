# Check whether quotations occur verbatim in a text after normalizing ligatures, quotation marks, hyphenation at line ends and white space.
import sys, re, unicodedata
def norm(s):
    s = s.replace("ﬁ", "fi").replace("ﬂ", "fl").replace("ﬀ", "ff").replace("ﬃ", "ffi").replace("ﬄ", "ffl")
    s = s.replace("“", '"').replace("”", '"').replace("‘", "'").replace("’", "'")
    s = s.replace("−", "-").replace("–", "-").replace("—", "-")
    s = re.sub(r"\s+", " ", s)
    return s.strip()
txt = norm(open(sys.argv[1], encoding="utf-8", errors="replace").read())
txt_nohyph = txt  # hyphenation at line ends cannot be undone reliably; report near matches instead
for q in open(sys.argv[2], encoding="utf-8").read().split("\n"):
    q = q.strip()
    if not q:
        continue
    nq = norm(q)
    i = txt.find(nq)
    if i >= 0:
        print("FOUND   :", q)
    else:
        # try ignoring spaces and hyphens entirely
        t2 = re.sub(r"[\s-]", "", txt); q2 = re.sub(r"[\s-]", "", nq)
        j = t2.find(q2)
        print("NOT VERBATIM" if j < 0 else "FOUND (ignoring spaces and hyphens):", q)
