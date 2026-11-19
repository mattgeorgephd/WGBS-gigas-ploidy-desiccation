# Fetch the full-text XML of an article from Europe PMC (open-access subset or
# author manuscripts where the service provides them) and write a plain-text
# version, one paragraph per line. Sends only the PMC identifier.
# The XML is untrusted data: it is parsed as text only.
# Usage: python3 -I epmc_fulltext.py <outdir> <label>=<PMCID> [...]
import html
import os
import re
import subprocess
import sys
import time

UA = "literature-metadata-check/1.0"


def fetch(url, dest, tries=3, pause=8):
    code = "000"
    for i in range(tries):
        r = subprocess.run(
            ["curl", "-sS", "-m", "90", "-A", UA, "-o", dest, "-w", "%{http_code}", url],
            capture_output=True, text=True)
        code = r.stdout.strip()
        if code == "200" or code in ("403", "404", "407"):
            return code
        time.sleep(pause * (i + 1))
    return code


def to_text(xml):
    xml = re.sub(r"(?is)<ref-list.*?</ref-list>", " ", xml)
    out = []
    # titles and paragraphs in document order
    for m in re.finditer(r"(?is)<(title|p|label)(?:\s[^>]*)?>(.*?)</\1>", xml):
        tag, body = m.group(1).lower(), m.group(2)
        body = re.sub(r"(?is)<xref[^>]*>(.*?)</xref>", r"\1", body)
        t = html.unescape(re.sub(r"(?s)<[^>]+>", "", body))
        t = re.sub(r"\s+", " ", t).strip()
        if not t:
            continue
        out.append(("## " + t) if tag == "title" else t)
    return "\n".join(out)


def main():
    outdir = sys.argv[1]
    os.makedirs(os.path.join(outdir, "epmc_fulltext"), exist_ok=True)
    textdir = os.path.join(os.path.dirname(outdir.rstrip("/")), "text")
    os.makedirs(textdir, exist_ok=True)
    for item in sys.argv[2:]:
        label, pmcid = item.split("=", 1)
        dest = os.path.join(outdir, "epmc_fulltext", label + "_" + pmcid + ".xml")
        code = fetch("https://www.ebi.ac.uk/europepmc/webservices/rest/" + pmcid + "/fullTextXML", dest)
        size = os.path.getsize(dest) if os.path.exists(dest) else 0
        print(label, pmcid, "http", code, "bytes", size)
        if code == "200" and size > 2000:
            raw = open(dest, encoding="utf-8", errors="replace").read()
            txt = to_text(raw)
            tp = os.path.join(textdir, label + "_" + pmcid + ".txt")
            open(tp, "w", encoding="utf-8").write(txt)
            print("   text:", tp, "lines", txt.count("\n") + 1)
        time.sleep(3)


if __name__ == "__main__":
    main()
