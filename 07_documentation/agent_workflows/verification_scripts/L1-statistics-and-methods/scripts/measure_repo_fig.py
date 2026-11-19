# Measure the bar heights of the simulation figure rendered in the authors' code repository (lossless PNG, 1920 x 1920).
# Same method as measure_fig2.py: gridlines give the axis (top = 1.00, bottom = 0.00); a bar's value is at the centre of its top outline stroke.
import sys
import numpy as np
from PIL import Image

im = np.asarray(Image.open(sys.argv[1]).convert("RGB")).astype(float)
H, W, _ = im.shape
bright = im.mean(axis=2)
cols = {"DSS": (230, 159, 0), "limma": (86, 180, 233), "methylKit Chisq (default)": (0, 158, 115), "methylKit Chisq OC": (240, 228, 66), "methylKit F OC": (0, 114, 178)}
cols = {k: np.array(v, dtype=float) for k, v in cols.items()}

def is_col(px, c, tol=30):
    return np.abs(px - c).sum(axis=-1) < tol

panels = {"a sensitivity": (117, 945, 50, 880), "b specificity": (1077, 1905, 50, 880), "c F-score": (110, 945, 1010, 1840)}
for name, (x0, x1, y0, y1) in panels.items():
    strip_x = slice(x1 - 30, x1 - 12)
    prof = bright[y0:y1, strip_x].mean(axis=1)
    grey = np.where((prof > 228) & (prof < 242))[0]
    gy0, gy1 = grey.min() + y0, grey.max() + y0
    rows = [y for y in range(gy0, gy1 + 1) if bright[y, strip_x].mean() > 247]
    groups = []
    for r in rows:
        if groups and r - groups[-1][-1] <= 1:
            groups[-1].append(r)
        else:
            groups.append([r])
    centres = []
    for g in groups:
        ys = np.arange(g[0] - 1, g[-1] + 2)
        w = np.array([bright[y, strip_x].mean() - 235.0 for y in ys]).clip(0)
        centres.append(float((w * ys).sum() / w.sum()))
    top, bottom = centres[0], centres[-1]
    print("\npanel %s: %d gridlines at %s; 1.00 at %.2f, 0.00 at %.2f; one pixel = %.5f" % (name, len(centres), [round(c, 1) for c in centres], top, bottom, 1.0 / (bottom - top)))
    sub = im[y0:y1, x0:x1]
    for tool, c in cols.items():
        mask = is_col(sub, c)
        colsum = mask.sum(axis=0)
        xs = np.where(colsum >= 2)[0]
        bars = []
        for x in xs:
            if bars and x - bars[-1][-1] <= 2:
                bars[-1].append(x)
            else:
                bars.append([x])
        bars = [b for b in bars if len(b) >= 6]
        vals = []
        for b in bars:
            xm = x0 + b[len(b)//2]
            ys = np.where(mask[:, b[len(b)//2 - 2: len(b)//2 + 3]].any(axis=1))[0]
            yfill = ys.min() + y0
            yy = np.arange(yfill - 7, yfill)
            dark = (200.0 - bright[yy, xm-2:xm+3].mean(axis=1)).clip(0)
            ytop = float((dark * yy).sum() / dark.sum()) if dark.sum() > 0 else float(yfill)
            vals.append((bottom - ytop) / (bottom - top))
        print("  %-26s %d bars: %s" % (tool, len(vals), ", ".join("%.4f" % v for v in vals)))
