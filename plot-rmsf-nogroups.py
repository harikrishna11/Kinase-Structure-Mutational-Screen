#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import sys

if len(sys.argv) < 4:
    print >>sys.stderr, "Usage: plot-rmsf.py in.xvg sample_name residue_offset [out.imgfile]"
    sys.exit(1)

vals = []
residue_offset = int(sys.argv[3])

with open(sys.argv[1], "r") as f:
    for l in f:
        if l.startswith("#") or l.startswith("@") or len(l.strip()) == 0:
            continue
        bits = l.split()
        vals.append(float(bits[1]))


fig, ax = plt.subplots(figsize = (50, 4))

group_offsets = numpy.arange(len(vals))
barwidth = 0.3

ax.bar(group_offsets, vals, barwidth, color = ["lightgray", "gray"], linewidth=0.5)

ax.set_xlabel("Residue")
ax.set_title("RMSF for sample " + sys.argv[2])
ax.set_ylabel("RMSF (nm)")
ax.set_xticks(group_offsets + (barwidth * 0.5))
ax.set_xticklabels(range(residue_offset, residue_offset + len(vals)), rotation = "vertical", size = "x-small")
ax.set_xlim(right=len(vals))
for tick in ax.get_xticklines():
    tick.set_visible(False)

if len(sys.argv) < 4:
    plt.show()
else:
    plt.savefig(sys.argv[4])
