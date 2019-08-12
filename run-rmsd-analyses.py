#!/usr/bin/env python

import sys
import read_groups
import subprocess
import copy
import os.path

if len(sys.argv) < 3:
    print >>sys.stderr, "Usage: run-rmsd-analyses.py outdir prefix analysiscmd analysiscmdarg ..."
    sys.exit(1)

outdir = sys.argv[1]
prefix = sys.argv[2]
cmd = sys.argv[3:]

# Pass dummy args, since I don't care about their offsets here:
groups = read_groups.get_rmsd_groups(*([1] * 11))

def make_proc(gname):

    filename = "%s_%s" % (prefix, gname)

    thiscmd = cmd + ["-o", os.path.join(outdir, "%s.xvg" % filename)]
    logout = os.path.join(outdir, "%s.out" % filename)
    logerr = os.path.join(outdir, "%s.err" % filename)
    with open(logout, "w") as fo, open(logerr, "w") as fe:
        proc = subprocess.Popen(thiscmd, stdout = fo, stderr = fe, stdin = subprocess.PIPE)
    proc.stdin.write("%s\n%s\n" % ("C-alpha", gname))
    proc.stdin.close()
    return proc

procs = map(make_proc, [g[0] for g in groups])

anyerr = False

for p in procs:
    if p.wait() != 0:
        anyerr = True

if anyerr:
    raise Exception("At least some proc failed!")
