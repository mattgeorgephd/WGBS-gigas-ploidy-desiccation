# Measure the bar heights of Figure 2 (panels a sensitivity, b specificity, c F-score) of
# Wreczycka et al. 2017 from the 300 dpi raster embedded in the publisher's PDF.
# The values are read from pixels, so they are approximate (about +-0.005); the paper prints no numbers.
# The image is untrusted data and is only decoded as pixels.
# Usage: python3 -I measure_fig2_bars.py <fig-000.jpg>
import sys
import numpy as np
from PIL import Image

im = np.asarray(Image.open(sys.argv[1]).convert("RGB")).astype(int)
H, W, _ = im.shape
print("image", W, "x", H)

# Okabe-Ito colours used by the figure (legend order)
tools = [("DSS", (230, 159, 0)), ("limma", (86, 180, 233)), ("methylKit Chisq (default)", (0, 158, 115)),
         ("methylKit Chisq OC", (240, 228, 66)), ("methylKit F OC", (0, 114, 178))]


def is_col(px, rgb, tol=38):
    return (np.abs(px - np.array(rgb)).sum(axis=-1) < tol)


def panel_axis(x0, x1, y0, y1):
    """Rows of the major white gridlines inside a panel, read in a column strip without bars."""
    strip = im[y0:y1, x0 + 4:x0 + 14, :]
    white = (strip.min(axis=2) > 248).all(axis=1)
    rows = np.where(white)[0] + y0
    # group consecutive rows
    groups = []
    for r in rows:
        if groups and r - groups[-1][-1] <= 1:
            groups[-1].append(r)
        else:
            groups.append([r])
    return [(g[0], g[-1], len(g)) for g in groups]


panels = {"a sensitivity": (112, 890, 60, 550), "b specificity": (1020, 1798, 60, 550), "c F-score": (105, 890, 695, 1186)}
for name, (x0, x1, y0, y1) in panels.items():
    gl = panel_axis(x0, x1, y0, y1)
    # nine gridlines per panel: 0.00, 0.125, ..., 1.00 (major and minor); the first and the last are 1.00 and 0.00
    major = gl
    centers = [0.5 * (g[0] + g[1]) for g in major]
    print("\npanel", name, "| gridline groups (first row, last row, thickness):", gl)
    if len(centers) < 2:
        print("  could not find gridlines"); continue
    # the major gridlines are at 0.00, 0.25, 0.50, 0.75, 1.00 from bottom to top
    centers = sorted(centers)
    top, bottom = centers[0], centers[-1]
    n_major = len(centers)
    print("  gridline centres:", [float(round(c, 1)) for c in centers], "->", n_major, "lines, 1.00 at the top one and 0.00 at the bottom one")
    assert n_major == 9, "expected nine gridlines"
    scale = (bottom - top) / 1.0
    sub = im[y0:y1, x0:x1, :]
    for tname, rgb in tools:
        mask = is_col(sub, rgb)
        colsum = mask.sum(axis=0)
        xs = np.where(colsum > 0)[0]
        bars = []
        for x in xs:
            if bars and x - bars[-1][-1] <= 2:
                bars[-1].append(x)
            else:
                bars.append([x])
        bars = [b for b in bars if len(b) >= 6]
        vals = []
        for b in bars:
            mid = b[len(b) // 2 - 2: len(b) // 2 + 3]
            tops = []
            for x in mid:
                rr = np.where(mask[:, x])[0]
                if len(rr):
                    tops.append(rr.min() + y0)
            if tops:
                ytop = float(np.median(tops)) - 1.5  # black outline of about 3 px above the fill
                vals.append(float(round((bottom - ytop) / scale, 3)))
        print("  %-28s bars found: %d  values by effect size 5,10,15,20,25: %s" % (tname, len(bars), vals))
