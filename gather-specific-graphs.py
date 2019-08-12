#!/usr/bin/env python

import sys
import glob
import os.path
import os
import subprocess

if len(sys.argv) != 3:
    print >>sys.stderr, "Usage: gather-specific-graphs.py protein_dir out_dir"

mydir = os.path.dirname(os.path.realpath(__file__))

protein_dir = sys.argv[1]
out_dir = sys.argv[2]
subdirs = filter(lambda d: os.path.isdir(os.path.join(protein_dir, d)), os.listdir(protein_dir))

if len(subdirs) == 0 or "wild" not in subdirs:
    print >>sys.stderr, "Need at least one subdirectory including one named 'wild'"
    sys.exit(1)

all_xvgs = set()

for d in subdirs:
    all_xvgs |= set(map(os.path.basename, glob.glob(os.path.join(protein_dir, d, "400ns_analyses", "*.xvg"))))
    
all_xvgs = list(all_xvgs)

rmsd_xvgs = filter(lambda d: d.find("_to_") == -1, all_xvgs)
dist_xvgs = filter(lambda d: d.find("_to_") != -1, all_xvgs)

def dirname_to_sname(d):
    if d == "wild":
        return "WT"
    else:
        return d

procs = []

def plot_xvgs(xvgs, plotname, ylabel):
    
    for x in xvgs:
        metric_name = x.replace(".xvg", "")
        args = [os.path.join(mydir, "plot-mavg-lines.py"), "%s %s" % (metric_name, plotname), "Time (ns)", ylabel]
        for d in subdirs:
            path = os.path.join(protein_dir, d, "400ns_analyses", x)
            if os.path.exists(path):
                args.extend([path, dirname_to_sname(d)])
        args.append(os.path.join(out_dir, "%s.pdf" % metric_name))
        outfile = os.path.join(out_dir, "%s.out" % metric_name)
        errfile = os.path.join(out_dir, "%s.out" % metric_name)
        with open(outfile, "w") as fo, open(errfile, "w") as fe:
            procs.append(subprocess.Popen(args, stdout = fo, stderr = fe))
    
plot_xvgs(rmsd_xvgs, "RMSD vs. time", "RMSD")
plot_xvgs(dist_xvgs, "inter-atomic distance vs. time", "Distance (nm)")

anyerrors = False

for p in procs:
    if p.wait() != 0:
        anyerrors = True

if anyerrors:
    raise Exception("At least some process failed")
