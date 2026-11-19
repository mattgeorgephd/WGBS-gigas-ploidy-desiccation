"""Independent rebuild of the per-CpG coverage and methylated-count matrices.
One coverage file is held in memory at a time. CpG sites are indexed from my own scan of the reference.
Outputs (in work/): cg_gpos.npy (genome coordinate of every reference CpG), seq table, cov.npy, nC.npy (uint16, CpG x library),
first_idx.npy (reference CpG index of the first line of each file), per-library summary."""
import sys, os, time, json
import numpy as np, pandas as pd

GENOME = "/home/user/gannet_cache/genome/cgigas_uk_roslin_v1_genomic-mito.fa"
COVDIR = "/home/user/gannet_cache/cov_files"
OUT = sys.argv[1]
lib = pd.read_csv("/home/user/WGBS-gigas-ploidy-desiccation/03_input_files/metadata/wgbs_libraries.csv").sort_values("legacy_sample_index")
t0 = time.time()

# ---- reference CpG index -------------------------------------------------------------------
names, lengths, cg_list = [], [], []
def finish(name, chunks):
    seq = np.frombuffer(b"".join(chunks), dtype=np.uint8)
    up = seq & 0xDF                                   # ASCII upper case
    pos = np.flatnonzero((up[:-1] == ord("C")) & (up[1:] == ord("G")))   # 0-based position of the C
    names.append(name); lengths.append(len(seq)); cg_list.append(pos.astype(np.int64))
name, chunks = None, []
with open(GENOME, "rb") as fh:
    for line in fh:
        if line.startswith(b">"):
            if name is not None: finish(name, chunks)
            name = line[1:].split()[0].decode(); chunks = []
        else:
            chunks.append(line.rstrip(b"\r\n"))
finish(name, chunks)
lengths = np.array(lengths, dtype=np.int64)
offsets = np.concatenate([[0], np.cumsum(lengths)[:-1]])
cg_gpos = np.concatenate([p + o for p, o in zip(cg_list, offsets)])
cg_seq = np.concatenate([np.full(len(p), i, dtype=np.int16) for i, p in enumerate(cg_list)])
n_cpg = len(cg_gpos)
print(f"[{time.time()-t0:6.0f}s] reference: {len(names)} sequences, {lengths.sum()} bp, {n_cpg} CG dinucleotides (case-insensitive, forward strand)", flush=True)
assert np.all(np.diff(cg_gpos) > 0)
np.save(os.path.join(OUT, "cg_gpos.npy"), cg_gpos); np.save(os.path.join(OUT, "cg_seq.npy"), cg_seq)
pd.DataFrame({"name": names, "length": lengths, "offset": offsets, "n_cg": [len(p) for p in cg_list]}).to_csv(os.path.join(OUT, "sequences.csv"), index=False)
del cg_list
seq_index = {n: i for i, n in enumerate(names)}

# ---- coverage files, one at a time ----------------------------------------------------------
n_lib = len(lib)
cov = np.lib.format.open_memmap(os.path.join(OUT, "cov.npy"), mode="w+", dtype=np.uint16, shape=(n_cpg, n_lib), fortran_order=True)
nC = np.lib.format.open_memmap(os.path.join(OUT, "nC.npy"), mode="w+", dtype=np.uint16, shape=(n_cpg, n_lib), fortran_order=True)
first_idx = np.zeros(n_lib, dtype=np.int64)
rows = []
for j, (lid, fn) in enumerate(zip(lib.library_id, lib.cov_file)):
    d = pd.read_csv(os.path.join(COVDIR, fn), sep="\t", header=None, usecols=[0, 1, 2, 4, 5],
                    names=["chr", "start", "end", "pct", "nC", "nT"],
                    dtype={"chr": "category", "start": np.int64, "end": np.int64, "nC": np.int64, "nT": np.int64})
    cats = d["chr"].cat.categories
    cat_to_seq = np.array([seq_index[c] for c in cats])
    sidx = cat_to_seq[d["chr"].cat.codes.to_numpy()]
    start = d["start"].to_numpy()
    gkey = offsets[sidx] + start                       # start is the 0-based position of the C
    idx = np.searchsorted(cg_gpos, gkey)
    idx_c = np.minimum(idx, n_cpg - 1)
    hit = cg_gpos[idx_c] == gkey
    n_not_cg = int((~hit).sum())
    n_end_bad = int((d["end"].to_numpy() - start != 2).sum())
    c = d["nC"].to_numpy(); t = d["nT"].to_numpy(); cv = c + t
    assert cv.max() < 65535 and cv.min() >= 1, (lid, cv.min(), cv.max())
    n_dup = len(idx) - len(np.unique(idx[hit]))
    cov[idx[hit], j] = cv[hit].astype(np.uint16); nC[idx[hit], j] = c[hit].astype(np.uint16)
    first_idx[j] = idx[0] if hit[0] else -1
    rows.append(dict(library_id=lid, n_rows=len(d), n_rows_not_at_reference_CG=n_not_cg, n_rows_end_minus_start_not_2=n_end_bad,
                     n_duplicate_positions=n_dup - n_not_cg if n_not_cg else n_dup, max_coverage=int(cv.max()), mean_coverage=float(cv.mean()),
                     median_coverage=float(np.median(cv)), n_cov_ge10=int((cv >= 10).sum()),
                     first_line=f"{d['chr'].iloc[0]}:{start[0]}", first_line_coverage=int(cv[0]),
                     chr_blocks_sorted=bool(np.all(np.diff(gkey[np.concatenate([[True], np.diff(sidx) == 0])] ) != 0))))
    print(f"[{time.time()-t0:6.0f}s] {lid:10s} rows={len(d)} not_at_CG={n_not_cg} dup={n_dup} maxcov={cv.max()} median={np.median(cv)}", flush=True)
    del d, sidx, start, gkey, idx, idx_c, hit, c, t, cv
cov.flush(); nC.flush()
np.save(os.path.join(OUT, "first_idx.npy"), first_idx)
pd.DataFrame(rows).to_csv(os.path.join(OUT, "library_summary.csv"), index=False)
print(f"[{time.time()-t0:6.0f}s] done", flush=True)
