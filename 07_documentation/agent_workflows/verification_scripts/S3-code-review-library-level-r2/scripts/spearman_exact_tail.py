# Exact lower tail of Spearman's statistic S = sum(d^2) for n untied pairs, up to S = CAP.
# Dynamic programming over the set of ranks already given to positions 1..pos (bit mask); a mask is dropped when even
# the best completion (sorted matching, which minimises a sum of squared distances on the line) exceeds CAP.
# Counts are exact 64-bit integers; the result for a small n is checked against full enumeration.
import numpy as np, math, sys, itertools, time
def tail_counts(n, cap):
    ranks = np.arange(1, n + 1)
    masks = np.array([0], dtype=np.int64); A = np.zeros((1, cap + 1), dtype=np.int64); A[0, 0] = 1
    biggest = 1
    for pos in range(1, n + 1):
        cand = []
        for r in range(1, n + 1):
            bit = 1 << (r - 1)
            rows = np.nonzero((masks & bit) == 0)[0]
            if len(rows): cand.append((r, rows, masks[rows] | bit))
        allnew = np.unique(np.concatenate([c[2] for c in cand]))
        bits = ((allnew[:, None] >> np.arange(n)) & 1).astype(np.int64)
        cost_in = (bits * (np.cumsum(bits, axis=1) - ranks) ** 2).sum(axis=1)
        cost_out = ((1 - bits) * (pos + np.cumsum(1 - bits, axis=1) - ranks) ** 2).sum(axis=1)
        allnew = allnew[cost_in + cost_out <= cap]
        del bits
        biggest = max(biggest, len(allnew))
        if len(allnew) * (cap + 1) * 8 > 1.2e9: raise MemoryError("layer %d would need %d rows" % (pos, len(allnew)))
        B = np.zeros((len(allnew), cap + 1), dtype=np.int64)
        for r, rows, nm in cand:
            sh = (pos - r) ** 2
            if sh > cap: continue
            idx = np.searchsorted(allnew, nm); idx[idx == len(allnew)] = 0
            ok = allnew[idx] == nm
            B[idx[ok], sh:] += A[rows[ok], :cap + 1 - sh]
        masks, A = allnew, B
    assert len(masks) == 1 and masks[0] == (1 << n) - 1
    return [int(x) for x in A[0]], biggest
# check against full enumeration for n = 8
n = 8; full = {}
for p in itertools.permutations(range(n)):
    s = sum((i - p[i]) ** 2 for i in range(n)); full[s] = full.get(s, 0) + 1
cnt, _ = tail_counts(8, 60)
assert all(cnt[q] == full.get(q, 0) for q in range(61)), "DP disagrees with enumeration for n = 8"
print("check for n = 8 against all 40320 permutations: identical counts up to S = 60")
t0 = time.time()
for n, qs in ((10, [30]), (12, [50]), (22, [396, 418, 430, 436, 486, 508])):
    cnt, biggest = tail_counts(n, max(qs))
    for q in qs:
        lower = sum(cnt[:q + 1])
        print("n = %d  S = %d  rho = %.6f  permutations with sum(d^2) <= S: %d  two-sided exact p = 2 * %d / %d! = %.6e" % (
            n, q, 1 - 6 * q / (n ** 3 - n), lower, lower, n, 2 * lower / math.factorial(n)))
    print("   largest number of masks kept in a layer: %d; seconds so far %.0f" % (biggest, time.time() - t0))
