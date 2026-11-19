import json, sys, os
d = sys.argv[1]
for fn in sorted(os.listdir(d)):
    if not fn.endswith(".json"):
        continue
    try:
        with open(os.path.join(d, fn), encoding="utf-8") as fh:
            m = json.load(fh)["message"]
    except Exception as e:
        print(fn, "ERROR", e); continue
    au = m.get("author", [])
    names = ["%s %s" % (a.get("family", a.get("name", "?")), a.get("given", "")) for a in au]
    def dp(k):
        v = m.get(k, {}).get("date-parts", [[None]])[0]
        return "-".join(str(x) for x in v)
    print("==", fn)
    print("  DOI:", m.get("DOI"))
    print("  type:", m.get("type"))
    print("  title:", "; ".join(m.get("title", [])))
    print("  authors (%d):" % len(names), "; ".join(names))
    print("  container:", "; ".join(m.get("container-title", [])), "| short:", "; ".join(m.get("short-container-title", [])))
    print("  volume:", m.get("volume"), "issue:", m.get("issue"), "page:", m.get("page"), "article-number:", m.get("article-number"))
    print("  published:", dp("published"), "print:", dp("published-print"), "online:", dp("published-online"), "issued:", dp("issued"))
    print("  publisher:", m.get("publisher"))
    rel = m.get("relation", {})
    if rel:
        print("  relation:", {k: [x.get("id") for x in v] for k, v in rel.items()})
    upd = m.get("update-to")
    if upd:
        print("  update-to:", upd)
    if m.get("updated-by"):
        print("  updated-by:", m.get("updated-by"))
