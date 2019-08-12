#!/usr/bin/env python

import numpy
import sys

if len(sys.argv) < 3:
    print >>sys.stderr, "Usage: find-stable-region.py in.xvg window-size"
    sys.exit(1)

window = float(sys.argv[2])
times = []
vals = []

with open(sys.argv[1], "r") as f:

    for l in f:

        if l.startswith("@") or l.startswith("#"):
            continue
        if len(l.strip()) == 0:
            continue
        bits = l.split()
        times.append(float(bits[0]))
        vals.append(float(bits[1]))

vals = numpy.array(vals)
vars = []

print >>sys.stderr, "Loaded", len(times), "points"

lim_idx = 0

for i in range(len(times)):

    if i % 1000 == 0:
        print >>sys.stderr, "Progress: time", times[i]

    start_time = times[i]
    time_lim = start_time + window
    while lim_idx < len(times) and times[lim_idx] < time_lim:
        lim_idx += 1

    window_vals = vals[i:lim_idx]
    vars.append((start_time, numpy.var(window_vals)))

    if lim_idx == len(times):
        break

(stable_region_time, stable_region_variance) = min(vars, key = lambda x : x[1])

print >>sys.stderr, "Most stable region starts at", stable_region_time, "with variance", stable_region_variance
print stable_region_time
