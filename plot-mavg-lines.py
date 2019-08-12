#!/usr/bin/env python

import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import sys

if len(sys.argv) < 4:
    print >>sys.stderr, "Usage: plot-rmsd.py title xlabel ylabel in.xvg sample_name [in2.xvg sample_name_2] ... [out.imgfile]"
    sys.exit(1)

title = sys.argv[1]
xlabel = sys.argv[2]
ylabel = sys.argv[3]

args = sys.argv[4:]
if len(args) % 2 == 0:
    outfile = None
else:
    outfile = args[-1]
    args = args[:-1]

moving_avg_window = 500
downsample_factor = 100

series_files = [x for (i, x) in enumerate(args) if i % 2 == 0]
series_names = [x for (i, x) in enumerate(args) if i % 2 == 1]

colours = ["red", "green", "blue", "orange", "purple", "yellow", "grey"]
colours.extend(["dark %s" % c for c in colours])

if len(colours) < len(series_names):
    print >>sys.stderr, "Define more colours!"
    sys.exit(1)

colours = colours[:len(series_names)]

def movingaverage(interval, window_size):
    window= numpy.ones(int(window_size))/float(window_size)
    return numpy.convolve(interval, window, 'same')

def rup(num, den):
    return int(math.ceil(float(num) / den))

for (fname, sname, col) in zip(series_files, series_names, colours):

    xpoints = []
    ypoints = []

    with open(fname, "r") as f:
        for l in f:
            if l.startswith("#") or l.startswith("@") or len(l.strip()) == 0:
                continue
            bits = l.split()
            xpoints.append(float(bits[0]))
            ypoints.append(float(bits[1]))

    avg_ypoints = movingaverage(ypoints, 500)
    xpoints = [x / 1000 for x in xpoints]

    if downsample_factor != 1:
        xpointssub = [x for (i, x) in enumerate(xpoints) if i % downsample_factor == (downsample_factor / 2)]
        ymaxes = [max(ypoints[i * downsample_factor : (i + 1) * downsample_factor]) for i in range(rup(len(ypoints), downsample_factor))]
        ymins = [min(ypoints[i * downsample_factor : (i + 1) * downsample_factor]) for i in range(rup(len(ypoints), downsample_factor))]

        # Translate to rectangles:
        rect_xs = []
        rect_ymaxes = []
        rect_ymins = []
        for i in range(rup(len(xpoints), downsample_factor)):
            rect_xs.append(xpoints[i * downsample_factor])
            try:
                rect_xs.append(xpoints[(i + 1) * downsample_factor])
            except IndexError:
                rect_xs.append(xpoints[-1])
            rect_ymaxes.append(ymaxes[i])
            rect_ymaxes.append(ymaxes[i])
            rect_ymins.append(ymins[i])
            rect_ymins.append(ymins[i])

        plt.fill_between(rect_xs, rect_ymins, rect_ymaxes, color = col, linewidth = 0.0, alpha = 0.25 if len(series_names) > 1 else 1.0)
    else:
        plt.plot(xpoints, ypoints, label = "_nolegend_", color = col, alpha = 0.5 if len(series_names) > 1 else 1.0)

    plt.plot(xpoints, avg_ypoints, label = sname, color = col)

plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.title(title)

lgd = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=5)

if outfile is None:
    plt.show()
else:
    plt.savefig(outfile, bbox_extra_artists = [lgd], bbox_inches = 'tight')
