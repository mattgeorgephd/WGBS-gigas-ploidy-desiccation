# Mean coverage of the 363 mitochondrial CpGs relative to each library's median coverage of all tested CpGs
# (normalized counts of united_cpg_legacy.csv.gz), and quartiles of the per-CpG cohort difference in percent methylation.
# Streams the gzipped table once; small memory. Usage: python3 -I mito_coverage_ratio.py <united_cpg_legacy.csv.gz> <mito_per_cpg.csv>
import csv, gzip, sys
from collections import Counter
hist = [Counter() for _ in range(22)]
msum = [0.0] * 22; mn = 0
with gzip.open(sys.argv[1], "rt") as f:
    r = csv.reader(f); h = next(r)
    idx = [h.index("coverage%d" % j) for j in range(1, 23)]
    n = 0
    for row in r:
        n += 1
        mito = row[0] == "NC_001276.1"
        if mito: mn += 1
        for j, i in enumerate(idx):
            v = float(row[i]); hist[j][v] += 1
            if mito: msum[j] += v
def median(c, n):
    k = sorted(c.items()); half = (n + 1) / 2.0; cum = 0; prev = None
    for v, m in k:
        if cum + m >= half:
            if n % 2 == 0 and cum + m == n / 2: continue
            return v
        cum += m
    return None
import statistics
meds = []
for j in range(22):
    vals = []
    for v, m in hist[j].items(): vals.extend([v] * m)
    meds.append(statistics.median(vals))
print("CpGs", n, "mitochondrial", mn)
ratios = [msum[j] / mn / meds[j] for j in range(22)]
for j in range(22):
    print("library %2d  median coverage of tested CpGs %6.1f  mean coverage of mito CpGs %7.1f  ratio %5.1f" % (j + 1, meds[j], msum[j] / mn, ratios[j]))
print("heat libraries (1-10): ratio %.1f to %.1f ; control libraries (11-22): ratio %.1f to %.1f" % (min(ratios[:10]), max(ratios[:10]), min(ratios[10:]), max(ratios[10:])))
print("ratio of cohort means of the ratio: %.2f" % (statistics.mean(ratios[10:]) / statistics.mean(ratios[:10])))
d = sorted(float(x["diff"]) for x in csv.DictReader(open(sys.argv[2])))
q = statistics.quantiles(d, n=4, method="inclusive")
print("per-CpG difference control minus heat, points: min %.2f q1 %.2f median %.2f q3 %.2f max %.2f" % (d[0], q[0], q[1], q[2], d[-1]))
