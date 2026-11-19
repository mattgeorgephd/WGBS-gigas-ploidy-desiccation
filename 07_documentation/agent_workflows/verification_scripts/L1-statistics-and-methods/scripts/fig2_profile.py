# Print the vertical colour profile through the middle of one bar, to see how the outline and the fill appear in the JPEG.
import sys
import numpy as np
from PIL import Image
im = np.asarray(Image.open(sys.argv[1]).convert("RGB")).astype(int)
x = int(sys.argv[2]); ya = int(sys.argv[3]); yb = int(sys.argv[4])
for y in range(ya, yb):
    print(y, im[y, x-2:x+3].mean(axis=0).round(0))
