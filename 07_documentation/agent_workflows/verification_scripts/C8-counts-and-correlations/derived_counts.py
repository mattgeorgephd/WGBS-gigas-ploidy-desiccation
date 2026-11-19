import numpy as np, os, sys
W = sys.argv[1]; TARGET = 12794074
cov = np.load(os.path.join(W, "cov.npy"), mmap_mode="r"); nC = np.load(os.path.join(W, "nC.npy"), mmap_mode="r")
n = cov.shape[0]
tot = np.zeros(n, dtype=np.int64); totC = np.zeros(n, dtype=np.int64); anyC = np.zeros(n, dtype=bool); anyT = np.zeros(n, dtype=bool); mx = np.zeros(n, dtype=np.int64)
tot_sets = {k: np.zeros(n, dtype=np.int64) for k in ("heat", "control")}
for j in range(22):
    c = np.asarray(cov[:, j]).astype(np.int64); m = np.asarray(nC[:, j]).astype(np.int64)
    tot += c; totC += m; anyC |= m > 0; anyT |= (c - m) > 0; mx = np.maximum(mx, c)
    tot_sets["heat" if j < 10 else "control"] += c
res = {}
res["union (>=1 read in >=1 library)"] = int((tot > 0).sum())
res["CpGs with >=1 methylated read in any library"] = int(anyC.sum())
res["CpGs with >=1 unmethylated read in any library"] = int(anyT.sum())
res["CpGs with both a methylated and an unmethylated read (pooled)"] = int((anyC & anyT).sum())
for s in (2, 3, 4, 5, 10, 11, 12, 20, 22, 44):
    res[f"CpGs with pooled coverage (sum over 22 libraries) >= {s}"] = int((tot >= s).sum())
for s in (2, 3, 5, 10, 12):
    res[f"CpGs with pooled coverage >= {s} in heat libraries or in control libraries"] = int(((tot_sets['heat'] >= s) | (tot_sets['control'] >= s)).sum())
    res[f"CpGs with pooled coverage >= {s} in heat libraries and in control libraries"] = int(((tot_sets['heat'] >= s) & (tot_sets['control'] >= s)).sum())
res["CpGs covered in >=1 heat and >=1 control library"] = int(((tot_sets['heat'] > 0) & (tot_sets['control'] > 0)).sum())
for k, v in res.items(): print(f"{v:10d}  {v - TARGET:+9d}  {k}")
# the pooled-coverage threshold whose count brackets the target
cs = np.bincount(np.minimum(tot, 200))[::-1].cumsum()[::-1]
below = [(s, int(cs[s])) for s in range(1, 60) if abs(int(cs[s]) - TARGET) < 150000]
print("pooled-coverage thresholds with counts within 150,000 of the target:", below)
