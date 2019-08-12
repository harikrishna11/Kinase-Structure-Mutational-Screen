#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
    print >>sys.stderr, "Usage: show-xvg.py in.xvg [out.imgfile]"
    sys.exit(1)

xpoints = []
ypoints = []

with open(sys.argv[1], "r") as f:
    for l in f:
        if l.startswith("#") or l.startswith("@") or len(l.strip()) == 0:
            continue
        bits = l.split()
        xpoints.append(float(bits[0]))
        ypoints.append(float(bits[1]))

plt.plot(xpoints, ypoints)
if len(sys.argv) < 3:
    plt.show()
else:
    plt.savefig(sys.argv[2])
