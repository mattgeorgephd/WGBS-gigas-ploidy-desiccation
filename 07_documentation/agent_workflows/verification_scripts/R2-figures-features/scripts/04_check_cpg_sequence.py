# For CpG positions of the pipeline tables (0-based start of the C), read the reference at start-1 .. start+2
# and report what bases are there. Streams the FASTA; low memory.
import sys, gzip, csv, random, collections
fa = "/home/user/gannet_cache/genome/cgigas_uk_roslin_v1_genomic-mito.fa"
t05 = "/home/user/WGBS-gigas-ploidy-desiccation/04_output/05_differential_methylation/tables/"
want = collections.defaultdict(set)
n_dml = 0
for f in ["dml_ploidy_legacy.csv", "dml_heat_legacy.csv", "dml_ploidy_corrected.csv", "dml_heat_corrected.csv"]:
    with open(t05 + f) as h:
        r = csv.DictReader(h)
        for row in r:
            want[row["chr"]].add(int(row["start"])); n_dml += 1
random.seed(1)
samp = 0
with gzip.open(t05 + "united_cpg_legacy.csv.gz", "rt") as h:
    r = csv.reader(h); hdr = next(r); ic = hdr.index("chr"); isx = hdr.index("start")
    for row in r:
        if random.random() < 0.02:
            want[row[ic]].add(int(row[isx])); samp += 1
print("DML rows read:", n_dml, " random analysed CpGs added:", samp, " unique positions:", sum(len(v) for v in want.values()), file=sys.stderr)
res = collections.Counter(); lens = {}
def flush(name, chunks):
    if name is None: return
    seq = "".join(chunks); lens[name] = len(seq)
    for p in want.get(name, ()):
        di = seq[p:p+2].upper()
        res[di] += 1
        if di != "CG" and res["_shown"] < 10:
            res["_shown"] += 1; print("not CG:", name, p, seq[max(0,p-3):p+5], file=sys.stderr)
name = None; chunks = []
with open(fa) as h:
    for line in h:
        if line[0] == ">":
            flush(name, chunks); name = line[1:].split()[0]; chunks = []
        elif name in want:
            chunks.append(line.rstrip("\n"))
        else:
            chunks.append(line.rstrip("\n")) if False else None
flush(name, chunks)
print({k: v for k, v in res.items() if k != "_shown"})
