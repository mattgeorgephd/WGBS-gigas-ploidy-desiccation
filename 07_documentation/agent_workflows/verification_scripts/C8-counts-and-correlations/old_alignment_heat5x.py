"""Heat-only correlation at 5x from the OLD alignment run (fastp-trim.20201202), using the 5x bedgraph files on gannet
(chr, start, end, percent methylation; one file per heat library; zr3534_1's file is truncated on gannet).
Streams one file at a time; nothing is stored except one float32 column per library."""
import sys, os, subprocess, time, io
import numpy as np, pandas as pd
W = sys.argv[1]
base = "https://gannet.fish.washington.edu/panopea/WGBS-gigas-ploidy-desiccation/bed_files/RONIT/"
cg_gpos = np.load(os.path.join(W, "cg_gpos.npy")); seqs = pd.read_csv(os.path.join(W, "sequences.csv"))
seq_index = {n: i for i, n in enumerate(seqs.name)}; offsets = seqs.offset.to_numpy()
n_cpg = len(cg_gpos); t0 = time.time()
pm = np.full((n_cpg, 10), np.nan, dtype=np.float32, order="F")
info = []
for k in range(1, 11):
    url = f"{base}zr3534_{k}_R1.fastp-trim.20201202_5x.bedgraph"
    for attempt in range(3):
        p = subprocess.Popen(["curl", "-sS", "--fail", "--http1.1", "--max-time", "1800", url], stdout=subprocess.PIPE)
        try:
            d = pd.read_csv(p.stdout, sep="\t", header=None, names=["chr", "start", "end", "pct"],
                            dtype={"chr": "category", "start": np.int64, "end": np.int64, "pct": np.float64}, on_bad_lines="skip")
        except Exception as ex:
            print("  read error:", ex, flush=True); d = None
        rc = p.wait()
        if d is not None and (rc == 0 or k == 1): break
        print(f"  retry {k} (curl exit {rc})", flush=True)
    d = d.dropna()
    cats = d["chr"].cat.categories
    known = np.array([c in seq_index for c in cats])
    cat_to_seq = np.array([seq_index.get(c, -1) for c in cats])
    sidx = cat_to_seq[d["chr"].cat.codes.to_numpy()]
    ok = sidx >= 0
    gkey = offsets[np.where(ok, sidx, 0)] + d["start"].to_numpy()
    idx = np.minimum(np.searchsorted(cg_gpos, gkey), n_cpg - 1)
    hit = ok & (cg_gpos[idx] == gkey)
    pm[idx[hit], k - 1] = d["pct"].to_numpy()[hit].astype(np.float32)
    info.append(dict(library=f"zr3534_{k}", rows=len(d), rows_at_reference_CG=int(hit.sum()), curl_exit=rc,
                     last_row=f"{d['chr'].iloc[-1]}:{int(d['start'].iloc[-1])}", pct_min=float(d.pct.min()), pct_max=float(d.pct.max())))
    print(f"[{time.time()-t0:5.0f}s] {info[-1]}", flush=True)
    del d
np.save(os.path.join(W, "old_align_heat_5x_pm.npy"), pm)

def pearson(m):
    m = m.astype(np.float64); mu = m.mean(0); z = m - mu; C = z.T @ z; dd = np.sqrt(np.diag(C)); return C / np.outer(dd, dd)
def rng(r):
    off = r[np.triu_indices(len(r), 1)]; return f"min={off.min():.4f} max={off.max():.4f} median={np.median(off):.4f} (pairs={len(off)})"
have = ~np.isnan(pm)
# 9 complete libraries (2..10)
rows9 = np.flatnonzero(have[:, 1:].all(1)); r9 = pearson(pm[rows9][:, 1:])
print(f"old alignment, libraries 2-10, CpGs at >=5x in all nine: {len(rows9)}; Pearson (raw percent) {rng(r9)}")
rows10 = np.flatnonzero(have.all(1)); r10 = pearson(pm[rows10])
print(f"old alignment, all ten libraries, CpGs at >=5x in all ten (limited by the truncated zr3534_1 file): {len(rows10)}; Pearson {rng(r10)}")
# the same two computations with the current alignment (my matrices), raw percent methylation, no first-line loss
cov = np.load(os.path.join(W, "cov.npy"), mmap_mode="r"); nC = np.load(os.path.join(W, "nC.npy"), mmap_mode="r")
cur = np.full((n_cpg, 10), np.nan, dtype=np.float32, order="F")
for j in range(10):
    cv = np.asarray(cov[:, j]).astype(np.float32); nc = np.asarray(nC[:, j]).astype(np.float32)
    okc = cv >= 5; cur[okc, j] = 100 * nc[okc] / cv[okc]
hc = ~np.isnan(cur)
c9 = np.flatnonzero(hc[:, 1:].all(1)); print(f"current alignment, libraries 2-10, >=5x in all nine: {len(c9)}; Pearson (raw percent) {rng(pearson(cur[c9][:, 1:]))}")
c10 = np.flatnonzero(hc.all(1)); print(f"current alignment, all ten, >=5x in all ten: {len(c10)}; Pearson (raw percent) {rng(pearson(cur[c10]))}")
both = np.flatnonzero(hc[:, 1:].all(1) & have[:, 1:].all(1))
print(f"CpGs at >=5x in all nine libraries in both alignments: {len(both)}; per-library correlation old vs current percent:",
      " ".join(f"{np.corrcoef(pm[both, j].astype(np.float64), cur[both, j].astype(np.float64))[0,1]:.4f}" for j in range(1, 10)))
print("counts per library old (rows at >=5x) vs current (>=5x):", [(int(have[:, j].sum()), int(hc[:, j].sum())) for j in range(10)])
