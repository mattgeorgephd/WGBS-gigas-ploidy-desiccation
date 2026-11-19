"""Independent recomputation of the library correlations and the PCA from the rebuilt matrices.
Legacy view: no upper coverage filter; the CpG on the first line of a library's file is not counted for that library.
Normalization as methylKit::normalizeCoverage (median), in two floating-point evaluation orders:
  'step04'   : round((new_cov * nC) / cov)         (01_pipeline/04_cpg_coverage.R)
  'methylKit': round(new_cov * (nC / cov))         (methylKit 1.28.0 source)
and 'raw' = no normalization (100 * nC / cov)."""
import sys, os, time
import numpy as np, pandas as pd
W = sys.argv[1]
which = sys.argv[2].split(",") if len(sys.argv) > 2 else None
t0 = time.time()
lib = pd.read_csv("/home/user/WGBS-gigas-ploidy-desiccation/03_input_files/metadata/wgbs_libraries.csv").sort_values("legacy_sample_index").reset_index(drop=True)
ids = np.array(lib.library_id); n_lib = len(ids)
cohort = np.array(lib.cohort)
cov = np.load(os.path.join(W, "cov.npy"), mmap_mode="r"); nC = np.load(os.path.join(W, "nC.npy"), mmap_mode="r")
first_idx = np.load(os.path.join(W, "first_idx.npy"))
n_cpg = cov.shape[0]

def col(j, legacy=True):
    cv = np.array(cov[:, j], dtype=np.int64); nc = np.array(nC[:, j], dtype=np.int64)
    if legacy: cv[first_idx[j]] = 0; nc[first_idx[j]] = 0
    return cv, nc

def build(libs, x, order, legacy=True):
    """rows covered >= x in every library of libs; percent methylation matrix (float64)."""
    in_all = np.ones(n_cpg, dtype=bool); med = []
    for j in libs:
        cv, _ = col(j, legacy); ok = cv >= x
        med.append(np.median(cv[ok])); in_all &= ok
    rows = np.flatnonzero(in_all); med = np.array(med); sc = med.max() / med
    pm = np.empty((len(rows), len(libs)), dtype=np.float64)
    ndiff = 0
    for k, j in enumerate(libs):
        cv, nc = col(j, legacy); cv = cv[rows].astype(np.float64); nc = nc[rows].astype(np.float64)
        if order == "raw":
            pm[:, k] = 100.0 * nc / cv; continue
        new_cov = np.rint(sc[k] * cv)
        if order == "step04":
            new_c = np.rint(new_cov * nc / cv); new_t = np.rint(new_cov * (cv - nc) / cv)
        else:
            fC = nc / cv; fT = (cv - nc) / cv
            new_c = np.rint(new_cov * fC); new_t = np.rint(new_cov * fT)
        pm[:, k] = 100.0 * new_c / (new_c + new_t)
    return rows, pm, med, sc

def pearson(pm, chunk=500_000):
    n, p = pm.shape; mu = np.zeros(p)
    for s in range(0, n, chunk): mu += pm[s:s+chunk].sum(0)
    mu /= n; C = np.zeros((p, p))
    for s in range(0, n, chunk):
        z = pm[s:s+chunk] - mu; C += z.T @ z
    d = np.sqrt(np.diag(C)); return C / np.outer(d, d)

def rank_inplace(pm):
    for k in range(pm.shape[1]):
        u, inv, cnt = np.unique(pm[:, k], return_inverse=True, return_counts=True)
        cs = np.cumsum(cnt); pm[:, k] = (cs - (cnt - 1) / 2.0)[inv]

def pca(pm, chunk=500_000):
    n, p = pm.shape; sds = np.empty(n)
    for s in range(0, n, chunk): sds[s:s+chunk] = pm[s:s+chunk].std(axis=1, ddof=1)
    cut = np.quantile(sds, 0.5); keep = np.flatnonzero(sds > cut)
    G = np.zeros((p, p))
    for s in range(0, len(keep), chunk):
        a = pm[keep[s:s+chunk]]; z = (a - a.mean(1, keepdims=True)) / a.std(axis=1, ddof=1, keepdims=True); G += z.T @ z
    # G = X X^T with X = libraries x CpGs, each CpG centred and scaled across libraries (prcomp(t(x), center, scale.))
    w, v = np.linalg.eigh(G); o = np.argsort(w)[::-1]; w = np.clip(w[o], 0, None); v = v[:, o]
    return 100 * w / w.sum(), v * np.sqrt(w), len(keep), cut

def summarize(r, libs, label):
    iu = np.triu_indices(len(libs), 1); off = r[iu]
    co = cohort[libs]; same = (co[:, None] == co[None, :])[iu]
    s = f"{label}: min={off.min():.4f} max={off.max():.4f} median={np.median(off):.4f}"
    if same.any() and (~same).any(): s += f" | within-cohort median={np.median(off[same]):.4f} between-cohort median={np.median(off[~same]):.4f}"
    print("   " + s, flush=True); return off

sets = {"all": np.arange(n_lib), "heat": np.flatnonzero(lib.treatment == "heat"), "control": np.flatnonzero(lib.treatment == "control"),
        "diploid": np.flatnonzero(lib.ploidy == "diploid"), "triploid": np.flatnonzero(lib.ploidy == "triploid")}
jobs = [("all", 1), ("all", 10), ("all", 3), ("all", 5), ("all", 20), ("heat", 1), ("heat", 5), ("heat", 10), ("control", 1), ("control", 5), ("control", 10),
        ("diploid", 1), ("triploid", 1), ("diploid", 10), ("triploid", 10)]
T = "/home/user/WGBS-gigas-ploidy-desiccation/04_output/04_cpg_coverage/tables/"
res = []
for sname, x in jobs:
    if which and f"{sname}{x}" not in which: continue
    libs = sets[sname]
    rows, pm, med, sc = build(libs, x, "step04")
    print(f"\n[{time.time()-t0:5.0f}s] set={sname} ({len(libs)} libraries) >= {x}x: {len(rows)} CpGs in all; medians={med.tolist()} NaN rows={int(np.isnan(pm).any(1).sum())}", flush=True)
    r = pearson(pm); off = summarize(r, libs, "Pearson, methylKit-style normalization (step 04 order)")
    rec = dict(sample_set=sname, min_coverage=x, n_cpg=len(rows), pearson_min=off.min(), pearson_max=off.max(), pearson_median=np.median(off))
    if sname == "all" and x in (1, 10):
        saved = pd.read_csv(T + f"library_correlation_{x}x_legacy.csv", index_col=0).values
        print(f"   max |r - saved pipeline matrix| = {np.abs(r - saved).max():.3e}")
        rec["max_abs_diff_vs_pipeline"] = float(np.abs(r - saved).max())
        ve, scores, nkeep, cut = pca(pm)
        pv = pd.read_csv(T + "pca_variance.csv"); pv = pv[(pv["mode"] == "legacy") & (pv.min_coverage == x)].sort_values("pc").percent_variance.values
        ps = pd.read_csv(T + f"pca_scores_{x}x_legacy.csv")
        d1 = min(np.abs(scores[:, 0] - ps.PC1.values).max(), np.abs(scores[:, 0] + ps.PC1.values).max())
        d2 = min(np.abs(scores[:, 1] - ps.PC2.values).max(), np.abs(scores[:, 1] + ps.PC2.values).max())
        print(f"   PCA: {nkeep} CpGs with SD above the median SD ({cut:.4f}); PC1 {ve[0]:.4f}% PC2 {ve[1]:.4f}% PC3 {ve[2]:.4f}% | pipeline PC1 {pv[0]:.4f}% PC2 {pv[1]:.4f}% | max |score diff| PC1 {d1:.3e} PC2 {d2:.3e}")
        rec.update(pc1=ve[0], pc2=ve[1], pipeline_pc1=pv[0], pipeline_pc2=pv[1])
        # cohort separation on the leading PCs
        for k in range(4):
            s = scores[:, k]; h = s[cohort == "heat_experiment"]; c = s[cohort == "control_experiment"]
            pooled = np.sqrt(((len(h)-1)*h.var(ddof=1) + (len(c)-1)*c.var(ddof=1)) / (len(h)+len(c)-2))
            print(f"      PC{k+1} ({ve[k]:.1f}%): heat mean {h.mean():8.1f} [{h.min():8.1f},{h.max():8.1f}]  control mean {c.mean():8.1f} [{c.min():8.1f},{c.max():8.1f}]  standardized difference {(h.mean()-c.mean())/pooled:5.2f}")
    rank_inplace(pm); rs = pearson(pm); offs = summarize(rs, libs, "Spearman, same matrix")
    rec.update(spearman_min=offs.min(), spearman_max=offs.max(), spearman_median=np.median(offs))
    del pm
    rows, pm, _, _ = build(libs, x, "methylKit"); r2 = pearson(pm); off2 = summarize(r2, libs, "Pearson, normalization in methylKit's evaluation order")
    print(f"   max |r(step04 order) - r(methylKit order)| = {np.abs(r - r2).max():.3e}")
    rec.update(pearson_mk_min=off2.min(), pearson_mk_max=off2.max(), max_abs_diff_orders=float(np.abs(r - r2).max()))
    del pm
    rows, pm, _, _ = build(libs, x, "raw"); r3 = pearson(pm); off3 = summarize(r3, libs, "Pearson, no normalization (100*nC/cov)")
    rec.update(pearson_raw_min=off3.min(), pearson_raw_max=off3.max())
    if sname == "all" and x == 1:
        np.save(os.path.join(W, "pm_raw_rows_1x.npy"), rows)
    del pm
    res.append(rec)
    pd.DataFrame(res).to_csv(os.path.join(W, f"my_correlations{'_' + '_'.join(which) if which else ''}.csv"), index=False)
print(f"\n[{time.time()-t0:5.0f}s] done")
