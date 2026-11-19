# Independent measurement of the bar heights of Figure 2 of Wreczycka et al. 2017 from the image embedded in the PDF.
# Reads one PNG (argument 1); prints, per panel and tool, the bar heights on the 0-1 axis.
# Axis: the horizontal gridlines of each panel (white on grey) are located in a strip without bars; the top one is 1.00
# and the bottom one 0.00 (nine lines, every 0.125). Bar: the top edge is the centre of the dark outline stroke above the fill.
import sys
import numpy as np
from PIL import Image

im = np.asarray(Image.open(sys.argv[1]).convert("RGB")).astype(float)
H, W, _ = im.shape
print("image", W, "x", H)
bright = im.mean(axis=2)

legend = {"DSS": (1240, 882), "limma": (1240, 926), "methylKit Chisq (default)": (1240, 970), "methylKit Chisq OC": (1240, 1013), "methylKit F OC": (1240, 1057)}
cols = {k: np.median(im[y-8:y+8, x-8:x+8].reshape(-1, 3), axis=0) for k, (x, y) in legend.items()}

def is_col(px, c, tol=45):
    return np.abs(px - c).sum(axis=-1) < tol

# panel boxes (x0, x1, y0, y1) chosen generously; the grey panel is found inside
panels = {"a sensitivity": (100, 900, 40, 560), "b specificity": (1005, 1800, 40, 560), "c F-score": (95, 900, 675, 1200)}
for name, (x0, x1, y0, y1) in panels.items():
    # grey panel background: value near 235. Find its vertical extent in a column strip at the right end of the panel
    strip_x = slice(x1 - 40, x1 - 20)
    prof = bright[y0:y1, strip_x].mean(axis=1)
    grey = np.where((prof > 228) & (prof < 242))[0]
    gy0, gy1 = grey.min() + y0, grey.max() + y0
    # gridlines: local maxima above 247 inside the grey panel
    rows = [y for y in range(gy0, gy1 + 1) if bright[y, strip_x].mean() > 247]
    groups = []
    for r in rows:
        if groups and r - groups[-1][-1] <= 1:
            groups[-1].append(r)
        else:
            groups.append([r])
    centres = []
    for g in groups:
        w = np.array([bright[y, strip_x].mean() - 235.0 for y in range(g[0] - 1, g[-1] + 2)])
        ys = np.arange(g[0] - 1, g[-1] + 2)
        centres.append(float((w.clip(0) * ys).sum() / w.clip(0).sum()))
    top, bottom = centres[0], centres[-1]
    print("\npanel %s: grey panel rows %d-%d; %d gridlines at %s" % (name, gy0, gy1, len(centres), [round(c, 1) for c in centres]))
    print("  1.00 at %.1f, 0.00 at %.1f; one pixel = %.4f" % (top, bottom, 1.0 / (bottom - top)))
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
            # darkness-weighted centre of the outline stroke in the 6 rows above the first fill row
            yy = np.arange(yfill - 6, yfill + 1)
            dark = (200.0 - bright[yy, xm-2:xm+3].mean(axis=1)).clip(0)
            # exclude the fill rows themselves (they are not darker than 200 minus fill brightness for bright fills, but may be for dark fills)
            dark[-1] = 0.0
            if dark.sum() > 0:
                ytop = float((dark * yy).sum() / dark.sum())
            else:
                ytop = float(yfill)
            vals.append((bottom - ytop) / (bottom - top))
        print("  %-26s %d bars: %s" % (tool, len(vals), ", ".join("%.3f" % v for v in vals)))
