"""Counts from the independently rebuilt matrices: union, per-set counts, k-of-n counts, and an exhaustive
search over every subset of the 22 libraries for a union equal to the manuscript's 12,794,074."""
import sys, os, time, json
import numpy as np, pandas as pd
W = sys.argv[1]
TARGET = 12794074
t0 = time.time()
lib = pd.read_csv("/home/user/WGBS-gigas-ploidy-desiccation/03_input_files/metadata/wgbs_libraries.csv").sort_values("legacy_sample_index").reset_index(drop=True)
ids = list(lib.library_id); n_lib = len(ids)
cov = np.load(os.path.join(W, "cov.npy"), mmap_mode="r")
first_idx = np.load(os.path.join(W, "first_idx.npy"))
cg_seq = np.load(os.path.join(W, "cg_seq.npy"))
seqs = pd.read_csv(os.path.join(W, "sequences.csv"))
n_cpg = cov.shape[0]
print("reference CpG sites:", n_cpg, " libraries:", n_lib, " first-line CpG found for all:", bool((first_idx >= 0).all()))
sets = {"all": np.arange(n_lib), "diploid": np.flatnonzero(lib.ploidy == "diploid"), "triploid": np.flatnonzero(lib.ploidy == "triploid"),
        "control": np.flatnonzero(lib.treatment == "control"), "heat": np.flatnonzero(lib.treatment == "heat")}
setmask = {k: int(sum(1 << int(j) for j in v)) for k, v in sets.items()}
FULL = (1 << n_lib) - 1

def masks(x, legacy):
    m = np.zeros(n_cpg, dtype=np.uint32)
    for j in range(n_lib):
        col = np.asarray(cov[:, j]) >= x
        if legacy: col[first_idx[j]] = False
        m |= (col.astype(np.uint32) << np.uint32(j))
    return m

pipe = pd.read_csv("/home/user/WGBS-gigas-ploidy-desiccation/04_output/04_cpg_coverage/tables/cpg_counts_by_coverage.csv")
pipe = pipe[pipe["mode"] == "legacy"].set_index(["sample_set", "min_coverage"])
out_rows = []; kofn_rows = []
natural = {}
# natural subsets for the union search
for j in range(n_lib): natural[FULL & ~(1 << j)] = f"all but {ids[j]}"
for a in range(n_lib):
    for b in range(a + 1, n_lib): natural[FULL & ~(1 << a) & ~(1 << b)] = f"all but {ids[a]},{ids[b]}"
for k in range(1, n_lib + 1): natural[(1 << k) - 1] = f"first {k} libraries in legacy order"
for k in range(1, n_lib): natural[FULL & ~((1 << k) - 1)] = f"last {n_lib-k} libraries in legacy order"
grp = {}
for p in ("diploid", "triploid"):
    for t in ("control", "heat"):
        grp[f"{p}/{t}"] = int(sum(1 << int(j) for j in np.flatnonzero((lib.ploidy == p) & (lib.treatment == t))))
gnames = list(grp)
for r in range(1, 16):
    mk = 0; nm = []
    for i, g in enumerate(gnames):
        if r >> i & 1: mk |= grp[g]; nm.append(g)
    natural.setdefault(mk, "groups " + "+".join(nm))
shock = int(sum(1 << int(j) for j in np.flatnonzero(lib.acute_heat_shock == "yes")))
natural.setdefault(FULL & ~shock, "all but the four acute-heat-shock libraries")
natural.setdefault(setmask["control"] | (setmask["heat"] & ~shock), "controls + heat without acute shock")

for legacy in (False, True):
    for x in (1, 2, 3, 5, 10, 20):
        if legacy and x == 2: continue
        m = masks(x, legacy)
        tag = f"{'legacy(first line lost)' if legacy else 'all lines'} >= {x}x"
        pc = np.zeros(n_cpg, dtype=np.uint8)
        for j in range(n_lib): pc += ((m >> np.uint32(j)) & np.uint32(1)).astype(np.uint8)
        hist = np.bincount(pc, minlength=n_lib + 1)
        atleast = np.cumsum(hist[::-1])[::-1]             # atleast[k] = CpGs covered in >= k libraries
        for k in range(1, n_lib + 1):
            kofn_rows.append(dict(view=tag, min_coverage=x, legacy=legacy, k=k, cpg_in_at_least_k_libraries=int(atleast[k]), diff_from_target=int(atleast[k]) - TARGET))
        for s, sm in setmask.items():
            sel = m & np.uint32(sm)
            n_all = int((sel == np.uint32(sm)).sum()); n_any = int((sel != 0).sum())
            row = dict(view=tag, min_coverage=x, legacy=legacy, sample_set=s, in_all=n_all, in_any=n_any)
            if legacy and (s, x) in pipe.index:
                row["pipeline_in_all"] = int(pipe.loc[(s, x), "cpg_in_all_libraries"]); row["pipeline_in_any"] = int(pipe.loc[(s, x), "cpg_in_any_library"])
                row["match"] = (row["pipeline_in_all"] == n_all) and (row["pipeline_in_any"] == n_any)
            out_rows.append(row)
        print(f"[{time.time()-t0:5.0f}s] {tag}: any={int(atleast[1])} all22={int(atleast[n_lib])}  >=2 libs={int(atleast[2])}  >=3 libs={int(atleast[3])}", flush=True)
        if x == 1 and not legacy:
            # by sequence class
            cls = np.where(seqs.name.str.startswith("NC_0475"), 0, np.where(seqs.name == "NC_001276.1", 2, 1))[cg_seq]
            for c, nm in enumerate(("linkage groups", "unplaced scaffolds", "mitochondrial genome")):
                print(f"      union at 1x on {nm}: {int(((m != 0) & (cls == c)).sum())} of {int((cls == c).sum())} reference CpGs")
            never = int((m == 0).sum()); print(f"      reference CpGs with no read in any library: {never}")
        if x in (1, 2, 3) and not legacy or (legacy and x == 1):
            # exhaustive subset search: union(S) for every subset S of libraries
            cnt = np.bincount(m, minlength=1 << n_lib).astype(np.int64)
            Z = cnt.copy()                                  # Z[S] = number of CpGs whose mask is a subset of S
            for b in range(n_lib):
                Zr = Z.reshape(-1, 2, 1 << b)
                Zr[:, 1, :] += Zr[:, 0, :]
            S_all = np.arange(1 << n_lib, dtype=np.int64)
            union = n_cpg - Z[FULL ^ S_all]
            assert union[FULL] == int(atleast[1]) and union[0] == 0
            hits = np.flatnonzero(union == TARGET)
            near = np.flatnonzero(np.abs(union - TARGET) <= 100)
            sizes = np.array([bin(int(h)).count("1") for h in hits])
            print(f"      subsets of the 22 libraries whose union equals {TARGET}: {len(hits)} (within +-100: {len(near)}); sizes of exact hits: {sorted(set(sizes.tolist()))}")
            nat_hits = [(natural[int(h)], int(union[h])) for h in hits if int(h) in natural]
            print(f"      of these, natural subsets (leave-one/two-out, prefixes, group unions, shock-excluded): {nat_hits}")
            # closest natural subsets
            nat = sorted(((abs(int(union[k]) - TARGET), int(union[k]), v) for k, v in natural.items()))[:6]
            print("      closest natural subsets:", [(v, u, u - TARGET) for d, u, v in nat])
            for h in hits[:5]:
                print("        e.g.", [ids[j] for j in range(n_lib) if int(h) >> j & 1])
            # how dense are union values around the target: number of distinct subsets per value in +-5000
            w = np.flatnonzero(np.abs(union - TARGET) <= 5000)
            print(f"      subsets with union within +-5000 of target: {len(w)} -> {len(w)/10001:.2f} subsets per integer value (chance hits expected)")
            lo = [(ids[j], int(union[FULL & ~(1 << j)])) for j in range(n_lib)]
            print("      leave-one-out unions:", lo)
            del cnt, Z, union
        del m, pc
pd.DataFrame(out_rows).to_csv(os.path.join(W, "my_counts_by_set.csv"), index=False)
pd.DataFrame(kofn_rows).to_csv(os.path.join(W, "my_counts_k_of_22.csv"), index=False)
df = pd.DataFrame(out_rows)
chk = df[df.legacy & df.pipeline_in_all.notna()]
print(f"\ncomparison with the pipeline table (legacy mode): {int(chk.match.sum())} of {len(chk)} (sample set x threshold) rows match on both counts")
print(chk[~chk.match.astype(bool)].to_string())
k = pd.DataFrame(kofn_rows)
print("\nclosest k-of-22 counts to the target:")
print(k.reindex(k.diff_from_target.abs().sort_values().index).head(8).to_string())
print("\nall-lines view, >=1x and >=2x, k-of-22:")
print(k[(~k.legacy) & (k.min_coverage.isin([1, 2]))].pivot(index="k", columns="min_coverage", values="cpg_in_at_least_k_libraries").to_string())
