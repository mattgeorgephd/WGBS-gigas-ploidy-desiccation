import json, subprocess, sys
def get(url):
    out = subprocess.run(["curl", "-sS", "-m", "40", url], capture_output=True, text=True).stdout
    return json.loads(out)
def listing(url):
    items = []
    while url:
        j = get(url); items += j["data"]; url = j["links"].get("next")
    return items
root = listing("https://api.osf.io/v2/nodes/ygctb/files/osfstorage/?page%5Bsize%5D=100")
for it in root:
    a = it["attributes"]; print("ROOT", a["kind"], a["name"])
for it in root:
    a = it["attributes"]
    if a["kind"] == "folder" and a["name"] in ("genome-feature-files", "code"):
        sub = listing(it["relationships"]["files"]["links"]["related"]["href"] + "?page%5Bsize%5D=100")
        for s in sub:
            b = s["attributes"]
            if a["name"] == "code" and "Feature" not in b["name"] and "feature" not in b["name"]: continue
            print(a["name"], b["kind"], b["name"], b.get("size"), (b.get("extra") or {}).get("hashes", {}).get("md5"), s["links"].get("download", ""))
