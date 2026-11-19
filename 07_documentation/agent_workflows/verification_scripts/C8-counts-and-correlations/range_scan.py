"""Scan natural library sets and coverage thresholds for a pairwise Pearson range that could be written as '0.87 to 0.93'.
Legacy view (no upper filter, first line of each file lost), methylKit-style median normalization within the set."""
import sys, os, time, math
import numpy as np, pandas as pd
W = sys.argv[1]; t0 = time.time()
lib = pd.read_csv("/home/user/WGBS-gigas-ploidy-desiccation/03_input_files/metadata/wgbs_libraries.csv").sort_values("legacy_sample_index").reset_index(drop=True)
ids = np.array(lib.library_id); n_lib = 22
cov = np.load(os.path.join(W, "cov.npy"), mmap_mode="r"); nC = np.load(os.path.join(W, "nC.npy"), mmap_mode="r"); first_idx = np.load(os.path.join(W, "first_idx.npy"))
n_cpg = cov.shape[0]
grp = np.array([f"{p}/{t}" for p, t in zip(lib.ploidy, lib.treatment)])
sets = {"all 22": np.arange(22), "heat (10)": np.flatnonzero(lib.treatment == "heat"), "control (12)": np.flatnonzero(lib.treatment == "control"),
        "diploid (11)": np.flatnonzero(lib.ploidy == "diploid"), "triploid (11)": np.flatnonzero(lib.ploidy == "triploid")}
for g in sorted(set(grp)): sets[f"{g} only"] = np.flatnonzero(grp == g)
sets["heat without acute shock (6)"] = np.flatnonzero((lib.treatment == "heat") & (lib.acute_heat_shock == "no"))
sets["controls + heat without acute shock (18)"] = np.flatnonzero((lib.treatment == "control") | (lib.acute_heat_shock == "no"))

def col(j):
    cv = np.array(cov[:, j], dtype=np.int64); nc = np.array(nC[:, j], dtype=np.int64); cv[first_idx[j]] = 0; nc[first_idx[j]] = 0; return cv, nc
def rmat(libs, x):
    in_all = np.ones(n_cpg, dtype=bool); med = []
    for j in libs:
        cv, _ = col(j); ok = cv >= x; med.append(np.median(cv[ok])); in_all &= ok
    rows = np.flatnonzero(in_all); med = np.array(med); sc = med.max() / med
    pm = np.empty((len(rows), len(libs)))
    for k, j in enumerate(libs):
        cv, nc = col(j); cv = cv[rows].astype(float); nc = nc[rows].astype(float)
        ncov = np.rint(sc[k] * cv); c = np.rint(ncov * nc / cv); t = np.rint(ncov * (cv - nc) / cv); pm[:, k] = 100 * c / (c + t)
    mu = pm.mean(0); C = np.zeros((len(libs),) * 2)
    for s in range(0, len(rows), 500_000):
        z = pm[s:s+500_000] - mu; C += z.T @ z
    d = np.sqrt(np.diag(C)); return len(rows), C / np.outer(d, d)
def fmt(off):
    lo, hi = off.min(), off.max()
    near = (round(lo + 1e-12, 2), round(hi + 1e-12, 2)); outw = (math.floor(lo * 100) / 100, math.ceil(hi * 100) / 100); trunc = (math.floor(lo * 100) / 100, math.floor(hi * 100) / 100)
    tag = []
    if near == (0.87, 0.93): tag.append("MATCH(nearest)")
    if outw == (0.87, 0.93): tag.append("MATCH(outward)")
    if trunc == (0.87, 0.93): tag.append("MATCH(truncated)")
    return lo, hi, near, outw, trunc, " ".join(tag)
out = []
for x in (1, 3, 5, 10, 20):
    for name, libs in sets.items():
        n, r = rmat(libs, x)
        iu = np.triu_indices(len(libs), 1); lo, hi, near, outw, trunc, tag = fmt(r[iu])
        out.append(dict(united_set=name, pairs="all pairs of the set", min_coverage=x, n_cpg=n, r_min=lo, r_max=hi, nearest=f"{near[0]:.2f}-{near[1]:.2f}", outward=f"{outw[0]:.2f}-{outw[1]:.2f}", flag=tag))
        print(f"[{time.time()-t0:5.0f}s] >= {x:2d}x {name:42s} n={n:8d} r {lo:.4f}..{hi:.4f}  nearest {near[0]:.2f}-{near[1]:.2f}  outward {outw[0]:.2f}-{outw[1]:.2f}  {tag}", flush=True)
        if name == "all 22":
            np.save(os.path.join(W, f"r_all22_{x}x.npy"), r)
            subs = {"heat pairs": lib.treatment == "heat", "control pairs": lib.treatment == "control", "diploid pairs": lib.ploidy == "diploid", "triploid pairs": lib.ploidy == "triploid"}
            for g in sorted(set(grp)): subs[f"{g} pairs"] = grp == g
            for sn, m in subs.items():
                k = np.flatnonzero(np.asarray(m)); rr = r[np.ix_(k, k)][np.triu_indices(len(k), 1)]; lo, hi, near, outw, trunc, tag = fmt(rr)
                out.append(dict(united_set="all 22", pairs=sn, min_coverage=x, n_cpg=n, r_min=lo, r_max=hi, nearest=f"{near[0]:.2f}-{near[1]:.2f}", outward=f"{outw[0]:.2f}-{outw[1]:.2f}", flag=tag))
                print(f"            within the 22-library set, {sn:24s} r {lo:.4f}..{hi:.4f}  nearest {near[0]:.2f}-{near[1]:.2f}  outward {outw[0]:.2f}-{outw[1]:.2f}  {tag}", flush=True)
pd.DataFrame(out).to_csv(os.path.join(W, "range_scan.csv"), index=False)
d = pd.DataFrame(out); print("\nrows scanned:", len(d), " with a flag:", int((d.flag != "").sum())); print(d[d.flag != ""].to_string())
