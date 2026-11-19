import json, re, glob
for p in ["code/4_WGBS_GOterm-Annotation-BOTH.ipynb", "code/3_WGBS_DML_feature_location-new.ipynb", "code/4_WGBS_GOterm-Annotation-RONIT.ipynb", "code/4_WGBS_GOterm-Annotation-YAAMINI.ipynb", "code/3_WGBS_DML_feature_location.ipynb"]:
    nb = json.load(open(p, encoding="utf-8")); cells = nb.get("cells", [])
    print("=====", p)
    for i, c in enumerate(cells):
        src = "".join(c.get("source", ""))
        if re.search(r"wc -l|union|unionbedg|1x", src):
            outs = []
            for o in c.get("outputs", []):
                if "text" in o: outs.append("".join(o["text"]))
                for k, v in (o.get("data") or {}).items():
                    if k.startswith("text"): outs.append("".join(v))
            out = "\n".join(outs)
            nums = [l.strip()[:200] for l in out.splitlines() if re.search(r"^\s*\d{4,}\b", l)]
            print(f"--- cell {i}: {src.strip()[:330].replace(chr(10), ' | ')}")
            for l in nums[:12]: print("      OUT:", l)
