import numpy as np, pandas as pd, itertools, sys
T = "/home/user/WGBS-gigas-ploidy-desiccation/04_output/04_cpg_coverage/tables/"
lib = pd.read_csv("/home/user/WGBS-gigas-ploidy-desiccation/03_input_files/metadata/wgbs_libraries.csv")
ids = list(lib.library_id)
cohort = np.array(lib.cohort); ploidy = np.array(lib.ploidy); shock = np.array(lib.acute_heat_shock)
group = np.array([f"{p}/{t}" for p, t in zip(lib.ploidy, lib.treatment)])
iu = np.triu_indices(22, 1)

def desc(v):
    return f"n={len(v):3d} min={v.min():.4f} max={v.max():.4f} median={np.median(v):.4f} mean={v.mean():.4f}"

for tag in ("1x_legacy", "1x_corrected", "10x_legacy", "10x_corrected"):
    m = pd.read_csv(T + f"library_correlation_{tag}.csv", index_col=0)
    assert list(m.index) == ids and list(m.columns) == ids
    r = m.values
    assert np.allclose(r, r.T) and np.allclose(np.diag(r), 1)
    off = r[iu]
    same_coh = (cohort[:, None] == cohort[None, :])[iu]
    print(f"\n######## {tag}")
    print("all pairs            ", desc(off))
    i, j = iu[0][off.argmin()], iu[1][off.argmin()]; print(f"   min pair: {ids[i]} - {ids[j]}")
    i, j = iu[0][off.argmax()], iu[1][off.argmax()]; print(f"   max pair: {ids[i]} - {ids[j]}")
    print("within cohort (both) ", desc(off[same_coh]))
    print("between cohorts      ", desc(off[~same_coh]))
    for c in ("heat_experiment", "control_experiment"):
        k = ((cohort[:, None] == c) & (cohort[None, :] == c))[iu]
        print(f"within {c:19s}", desc(off[k]))
    for p in ("diploid", "triploid"):
        k = ((ploidy[:, None] == p) & (ploidy[None, :] == p))[iu]
        print(f"within {p:19s}", desc(off[k]))
    k = (ploidy[:, None] != ploidy[None, :])[iu]
    print("between ploidies     ", desc(off[k]))
    for g in sorted(set(group)):
        k = ((group[:, None] == g) & (group[None, :] == g))[iu]
        print(f"within {g:19s}", desc(off[k]))
    k = (group[:, None] == group[None, :])[iu]
    print("within the 4 groups  ", desc(off[k]))
    print("between the 4 groups ", desc(off[~k]))
    # per-library mean correlation with others
    mean_other = (r.sum(1) - 1) / 21
    print("per-library mean r with the other 21:", " ".join(f"{i}={v:.3f}" for i, v in zip(ids, mean_other)))
    # share of pairs inside the manuscript's quoted range
    print(f"pairs with 0.865 <= r < 0.935 (would round into 0.87..0.93): {int(((off>=0.865)&(off<0.935)).sum())} of {len(off)};"
          f" below 0.865: {int((off<0.865).sum())}; at or above 0.935: {int((off>=0.935).sum())}")
    # largest subset of libraries whose pairwise r all round into [0.87, 0.93]: greedy + exhaustive over removal order
    ok = (r >= 0.865) & (r < 0.935)
    np.fill_diagonal(ok, True)
    best = []
    # exact maximum clique by simple branch and bound (22 nodes)
    def expand(cur, cand):
        global best
        if len(cur) > len(best): best = cur[:]
        if len(cur) + len(cand) <= len(best): return
        for idx, v in enumerate(cand):
            expand(cur + [v], [u for u in cand[idx+1:] if ok[v, u]])
    expand([], list(range(22)))
    sub = r[np.ix_(best, best)][np.triu_indices(len(best), 1)]
    print(f"largest set of libraries with every pairwise r in [0.865, 0.935): {len(best)} libraries: {[ids[b] for b in best]}")
    if len(best) > 1: print("   its range:", desc(sub))
    # subsets (maximal by size) whose min rounds to 0.87 and max rounds to 0.93 exactly
    found = None
    for size in range(22, 2, -1):
        if found: break
        if size > len(best): continue
        # search cliques of this size whose min in [0.865,0.875) and max in [0.925,0.935)
        cnt = 0
        def rec(cur, cand):
            global found, cnt
            if found: return
            if len(cur) == size:
                s = r[np.ix_(cur, cur)][np.triu_indices(size, 1)]
                if 0.865 <= s.min() < 0.875 and 0.925 <= s.max() < 0.935: found = cur[:]
                return
            if len(cur) + len(cand) < size: return
            for idx, v in enumerate(cand):
                rec(cur + [v], [u for u in cand[idx+1:] if ok[v, u]])
        rec([], list(range(22)))
    if found:
        s = r[np.ix_(found, found)][np.triu_indices(len(found), 1)]
        print(f"largest subset whose range rounds to exactly 0.87..0.93: {len(found)} libraries {[ids[b] for b in found]}: {desc(s)}")
    else:
        print("no subset of 3 or more libraries has a pairwise range that rounds to exactly 0.87..0.93")
