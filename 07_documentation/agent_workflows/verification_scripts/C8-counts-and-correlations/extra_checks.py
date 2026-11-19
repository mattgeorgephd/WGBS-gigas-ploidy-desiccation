import sys, os, numpy as np, pandas as pd
W = sys.argv[1]
GENOME = "/home/user/gannet_cache/genome/cgigas_uk_roslin_v1_genomic-mito.fa"
T = "/home/user/WGBS-gigas-ploidy-desiccation/04_output/04_cpg_coverage/tables/"
lib = pd.read_csv("/home/user/WGBS-gigas-ploidy-desiccation/03_input_files/metadata/wgbs_libraries.csv").sort_values("legacy_sample_index").reset_index(drop=True)
ids = list(lib.library_id); n_lib = 22; cohort = np.array(lib.cohort)
cov = np.load(os.path.join(W, "cov.npy"), mmap_mode="r")
n_cpg = cov.shape[0]

# 1. soft-mask class of every reference CpG, in reference order
cls_list = []
def finish(chunks):
    seq = np.frombuffer(b"".join(chunks), dtype=np.uint8); up = seq & 0xDF
    pos = np.flatnonzero((up[:-1] == ord("C")) & (up[1:] == ord("G")))
    c_low = seq[pos] >= 97; g_low = seq[pos + 1] >= 97
    cls_list.append(np.where(c_low & g_low, 1, np.where(c_low | g_low, 2, 0)).astype(np.uint8))   # 0 upper, 1 lower, 2 mixed
name, chunks = None, []
with open(GENOME, "rb") as fh:
    for line in fh:
        if line.startswith(b">"):
            if name is not None: finish(chunks)
            name = line; chunks = []
        else: chunks.append(line.rstrip(b"\r\n"))
finish(chunks)
cls = np.concatenate(cls_list); assert len(cls) == n_cpg
anycov = np.zeros(n_cpg, dtype=bool); nlib = np.zeros(n_cpg, dtype=np.uint8); all10 = np.ones(n_cpg, dtype=bool)
cut1 = []; cut2 = []
for j in range(n_lib):
    c = np.asarray(cov[:, j]); anycov |= c > 0; nlib += (c > 0); all10 &= c >= 10
    v = c[c >= 1].astype(np.float64); cut1.append(np.quantile(v, 0.999)); cut2.append(np.quantile(v[v >= 2], 0.999))
print("1. Reference CpGs by soft-mask class, and how many have a read in at least one library:")
for k, nm in enumerate(("both bases upper case", "both bases lower case (soft-masked)", "mixed case")):
    m = cls == k
    print(f"   {nm:38s} n={int(m.sum()):9d}  covered in >=1 library: {int((m & anycov).sum()):9d} ({100*(m & anycov).sum()/m.sum():.2f}%)  in all 22 at >=1x: {int((m & (nlib == 22)).sum()):8d}  in all 22 at >=10x: {int((m & all10).sum()):7d}")
print(f"   covered CpGs that are not fully upper case: {int((anycov & (cls != 0)).sum())} of {int(anycov.sum())} ({100*(anycov & (cls != 0)).sum()/anycov.sum():.1f}%)")

# 4. upper cutoffs: percentile base rows >= 1 read vs rows >= 2 reads
pc = pd.read_csv(T + "library_coverage.csv")
print("\n4. 99.9th percentile of coverage per library:")
print("   base rows >= 2 reads (step 04/05 'corrected'):", [round(float(x), 2) for x in cut2])
print("   pipeline library_coverage.csv upper_cutoff   :", [round(float(x), 2) for x in pc.upper_cutoff])
print("   base rows >= 1 read (methRead mincov = 1)    :", [round(float(x), 2) for x in cut1])
print("   max |mine(>=2) - pipeline| =", float(np.abs(np.array(cut2) - pc.upper_cutoff.values).max()), " libraries where the two bases give different cutoffs:", int((np.array(cut1) != np.array(cut2)).sum()))

# 5. pipeline in-all counts vs the legacy workbook
import openpyxl
wb = openpyxl.load_workbook("/home/user/WGBS-gigas-ploidy-desiccation/08_archive/analyses/coverage_summary.xlsx", data_only=True)
ws = wb["total_cpg"]; old = {(r[1].value, r[0].value): r[2].value for r in ws.iter_rows(min_row=2) if r[0].value is not None}
cc = pd.read_csv(T + "cpg_counts_by_coverage.csv"); cc = cc[cc["mode"] == "legacy"]
nm = {"all": "both"}; ok = 0
for _, r in cc.iterrows():
    key = (nm.get(r.sample_set, r.sample_set), r.min_coverage); ok += int(old[key] == r.cpg_in_all_libraries)
print(f"\n5. pipeline legacy in-all counts equal to analyses/coverage_summary.xlsx sheet total_cpg: {ok} of {len(cc)}")

# 2, 3, 6: correlation structure
rng = np.random.default_rng(1)
iu = np.triu_indices(22, 1)
depth = pc.set_index("library_id").loc[ids, "mean_coverage"].values
for tag in ("1x_legacy", "10x_legacy"):
    r = pd.read_csv(T + f"library_correlation_{tag}.csv", index_col=0).values; off = r[iu]
    same = (cohort[:, None] == cohort[None, :])[iu]
    obs_med = np.median(off[same]) - np.median(off[~same]); obs_mean = off[same].mean() - off[~same].mean()
    cnt_med = cnt_mean = 0; B = 100000
    for b in range(B):
        p = rng.permutation(cohort); s = (p[:, None] == p[None, :])[iu]
        cnt_med += (np.median(off[s]) - np.median(off[~s])) >= obs_med; cnt_mean += (off[s].mean() - off[~s].mean()) >= obs_mean
    hh = ((cohort[:, None] == "heat_experiment") & (cohort[None, :] == "heat_experiment"))[iu]
    cc_ = ((cohort[:, None] == "control_experiment") & (cohort[None, :] == "control_experiment"))[iu]
    print(f"\n2. {tag}: within-cohort minus between-cohort: median diff {obs_med:.4f} (permutation p = {(cnt_med+1)/(B+1):.4f}); mean diff {obs_mean:.4f} (p = {(cnt_mean+1)/(B+1):.4f}); "
          f"heat-heat median {np.median(off[hh]):.4f}, control-control {np.median(off[cc_]):.4f}, heat-control {np.median(off[~same]):.4f}")
    # 3. how much of the pairwise r is explained by library depth
    X = np.column_stack([np.ones(len(off)), np.log(depth[iu[0]]) + np.log(depth[iu[1]])])
    beta, res, *_ = np.linalg.lstsq(X, off, rcond=None); fit = X @ beta
    r2 = 1 - ((off - fit) ** 2).sum() / ((off - off.mean()) ** 2).sum()
    mean_other = (r.sum(1) - 1) / 21
    print(f"3. {tag}: R^2 of pairwise r on log depth_i + log depth_j = {r2:.3f}; correlation of per-library mean r with library mean coverage = {np.corrcoef(mean_other, depth)[0,1]:.3f}"
          f" (excluding zr3644_1 and zr3644_14: {np.corrcoef(np.delete(mean_other, [10, 17]), np.delete(depth, [10, 17]))[0,1]:.3f})")
for x in (1, 10):
    ps = pd.read_csv(T + f"pca_scores_{x}x_legacy.csv"); pv = pd.read_csv(T + "pca_variance.csv"); pv = pv[(pv["mode"] == "legacy") & (pv.min_coverage == x)]
    print(f"\n6. pipeline PCA scores, {x}x legacy:")
    for k in (1, 2, 3, 4):
        s = ps[f"PC{k}"].values; h = np.sort(s[ps.treatment == "heat"]); c = np.sort(s[ps.treatment == "control"])
        # number of controls outside the heat range, and best single-threshold classification accuracy
        best = 0
        for thr in np.sort(s):
            for sign in (1, -1):
                pred = (sign * s >= sign * thr); acc = ((pred == (ps.treatment == "heat").values).sum()); best = max(best, acc)
        d = np.sort(s[ps.ploidy == "diploid"]); t = np.sort(s[ps.ploidy == "triploid"]); bestp = 0
        for thr in np.sort(s):
            for sign in (1, -1):
                pred = (sign * s >= sign * thr); bestp = max(bestp, (pred == (ps.ploidy == "diploid").values).sum())
        print(f"   PC{k} ({pv[pv.pc == k].percent_variance.values[0]:.2f}%): heat [{h.min():.1f}, {h.max():.1f}]; controls {np.round(c, 1).tolist()}; controls inside the heat range: {int(((c >= h.min()) & (c <= h.max())).sum())}; best one-threshold split by treatment {best}/22, by ploidy {bestp}/22")
