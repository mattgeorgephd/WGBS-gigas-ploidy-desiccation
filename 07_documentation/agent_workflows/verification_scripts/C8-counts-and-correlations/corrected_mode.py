"""Corrected mode, recomputed independently: every line read; CpGs at or above the library's 99.9th percentile removed.
Two bases for the percentile: rows with >= 2 reads (what step 04 uses) and rows with >= 1 read (what the legacy coverage script imports)."""
import sys, os, numpy as np, pandas as pd
W = sys.argv[1]; T = "/home/user/WGBS-gigas-ploidy-desiccation/04_output/04_cpg_coverage/tables/"
cov = np.load(os.path.join(W, "cov.npy"), mmap_mode="r"); nC = np.load(os.path.join(W, "nC.npy"), mmap_mode="r"); n_cpg = cov.shape[0]
cut = {1: [], 2: []}
for j in range(22):
    c = np.asarray(cov[:, j]); v = c[c >= 1].astype(float)
    cut[1].append(np.quantile(v, 0.999)); cut[2].append(np.quantile(v[v >= 2], 0.999))
pipe = pd.read_csv(T + "cpg_counts_by_coverage.csv"); pipe = pipe[(pipe["mode"] == "corrected") & (pipe.sample_set == "all")].set_index("min_coverage")
def run(x, base, want_r):
    in_all = np.ones(n_cpg, dtype=bool); in_any = np.zeros(n_cpg, dtype=bool); med = []
    for j in range(22):
        c = np.asarray(cov[:, j]).astype(np.int64); ok = (c >= x) & (c < cut[base][j]); med.append(np.median(c[ok])); in_all &= ok; in_any |= ok
    rows = np.flatnonzero(in_all); out = dict(min_coverage=x, percentile_base=f">= {base} read(s)", in_all=len(rows), in_any=int(in_any.sum()))
    if want_r:
        med = np.array(med); sc = med.max() / med; pm = np.empty((len(rows), 22))
        for j in range(22):
            cv = np.asarray(cov[:, j])[rows].astype(float); nc = np.asarray(nC[:, j])[rows].astype(float)
            ncov = np.rint(sc[j] * cv); a = np.rint(ncov * nc / cv); b = np.rint(ncov * (cv - nc) / cv); pm[:, j] = 100 * a / (a + b)
        mu = pm.mean(0); C = np.zeros((22, 22))
        for s in range(0, len(rows), 500_000):
            z = pm[s:s+500_000] - mu; C += z.T @ z
        d = np.sqrt(np.diag(C)); r = C / np.outer(d, d); off = r[np.triu_indices(22, 1)]
        out.update(r_min=off.min(), r_max=off.max())
        if base == 2:
            saved = pd.read_csv(T + f"library_correlation_{x}x_corrected.csv", index_col=0).values; out["max_abs_diff_vs_pipeline_matrix"] = float(np.abs(r - saved).max())
    return out
res = []
for x in (1, 3, 5, 10, 20):
    for base in (2, 1):
        o = run(x, base, want_r=(x in (1, 10)))
        if base == 2: o["pipeline_in_all"] = int(pipe.loc[x, "cpg_in_all_libraries"]); o["pipeline_in_any"] = int(pipe.loc[x, "cpg_in_any_library"])
        res.append(o); print(o, flush=True)
pd.DataFrame(res).to_csv(os.path.join(W, "corrected_mode_check.csv"), index=False)
