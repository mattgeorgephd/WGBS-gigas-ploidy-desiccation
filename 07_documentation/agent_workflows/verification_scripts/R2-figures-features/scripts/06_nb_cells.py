import json, sys
nb = json.load(open(sys.argv[1]))
for i, c in enumerate(nb["cells"]):
    src = "".join(c["source"])
    if c["cell_type"] == "markdown":
        print(f"--- [{i}] MD: {src[:600]}")
    else:
        outs = []
        for o in c.get("outputs", []):
            if "text" in o: outs.append("".join(o["text"])[:400])
            elif "data" in o and "text/plain" in o["data"]: outs.append("".join(o["data"]["text/plain"])[:400])
        print(f"--- [{i}] CODE: {src[:900]}")
        if outs: print("    OUT: " + " | ".join(outs)[:500].replace("\n", "\n         "))
