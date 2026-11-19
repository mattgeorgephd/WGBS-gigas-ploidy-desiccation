# Measure the bars of the counts figure of the authors' repository (TP, FN, FP, TN per effect size and tool; lossless PNG).
# Axis: horizontal gridlines at 0, 500, ..., 2500 (major) located in a strip without bars inside the first facet.
import sys
import numpy as np
from PIL import Image
im = np.asarray(Image.open(sys.argv[1]).convert("RGB")).astype(float)
H, W, _ = im.shape
bright = im.mean(axis=2)
cols = {"DSS": (230, 159, 0), "limma": (86, 180, 233), "default": (0, 158, 115), "ChisqOC": (240, 228, 66), "FOC": (0, 114, 178)}
cols = {k: np.array(v, dtype=float) for k, v in cols.items()}
# facets: find grey panel column ranges at a row inside the plot area without bars (y = 80 is inside the panels, above 2500)
row = bright[22]            # a row inside the facet header strips (darker grey), above the facet labels
hdr = (row > 205) & (row < 225)
xs = np.where(hdr)[0]
facets = []
for x in xs:
    if facets and x - facets[-1][-1] <= 3:
        facets[-1].append(x)
    else:
        facets.append([x])
facets = [(f[0], f[-1]) for f in facets if len(f) > 100]
print("facets (x ranges):", facets)
# gridlines in a strip at the far left of the first facet (x0+2 .. x0+10): no bars there
x0, x1 = facets[0]
strip = slice(x0 + 2, x0 + 12)
rows = [y for y in range(60, 1840) if bright[y, strip].mean() > 247]
groups = []
for r in rows:
    if groups and r - groups[-1][-1] <= 1:
        groups[-1].append(r)
    else:
        groups.append([r])
cent = []
for g in groups:
    ys = np.arange(g[0] - 1, g[-1] + 2)
    w = np.array([bright[y, strip].mean() - 235.0 for y in ys]).clip(0)
    cent.append((float((w * ys).sum() / w.sum()), len(g)))
print("gridlines (centre, width):", [(round(c, 1), w) for c, w in cent])
major = [c for c, w in cent if w >= 2]
print("major gridlines:", [round(c, 1) for c in major])
ytop, ybot = major[0], major[-1]   # 2500 and 0
scale = 2500.0 / (ybot - ytop)
print("2500 at %.2f, 0 at %.2f; one pixel = %.3f counts" % (ytop, ybot, scale))
for fi, (fx0, fx1) in enumerate(facets):
    sub = im[:, fx0:fx1 + 1]
    out = {}
    for tool, c in cols.items():
        mask = np.abs(sub - c).sum(axis=-1) < 30
        colsum = mask[60:1840].sum(axis=0)
        xs = np.where(colsum >= 2)[0]
        bars = []
        for x in xs:
            if bars and x - bars[-1][-1] <= 1:
                bars[-1].append(x)
            else:
                bars.append([x])
        bars = [b for b in bars if len(b) >= 3]
        vals = []
        for b in bars:
            xm = fx0 + b[len(b)//2]
            ys = np.where(mask[60:1840, b[len(b)//2]])[0]
            yfill = ys.min() + 60
            yy = np.arange(yfill - 7, yfill)
            dark = (200.0 - bright[yy, xm]).clip(0)
            ytopbar = float((dark * yy).sum() / dark.sum()) if dark.sum() > 0 else float(yfill)
            vals.append(((fx0 + b[0]), round((ybot - ytopbar) * scale, 1)))
        out[tool] = vals
    print("facet", fi + 1, {k: [v[1] for v in vs] for k, vs in out.items()})
    print("        x positions", {k: [v[0] for v in vs] for k, vs in out.items()})
