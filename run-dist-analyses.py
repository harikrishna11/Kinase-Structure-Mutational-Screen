#!/usr/bin/env python

import sys
import read_groups
import subprocess
import copy
import os.path

if len(sys.argv) < 3:
    print >>sys.stderr, "Usage: run-rmsd-analyses.py outdir analysiscmd analysiscmdarg ..."
    sys.exit(1)

outdir = sys.argv[1]
cmd = sys.argv[2:]

compare_groups = [("VAIK_K_NZ", "bridge_E_CD"),
                  ("VAIK_K_NZ", "DFG_D_CG"),
                  ("DFG_F_O", "DFG_p2_N"),
                  ("HRD_R_NE", "DFG_p1_O"),
                  ("GxGxxG_5", "RSA"),
                  ("RS1", "RS2"),
                  ("RS2", "RS3"),
                  ("RS3", "RS4"),
                  ("RS4", "RSA")]

def make_proc(gpair):

    thiscmd = cmd + ["-select", 'com of group "%s" plus com of group "%s"' % gpair, "-oav", os.path.join(outdir, "%s_to_%s.xvg" % gpair)]
    logout = os.path.join(outdir, "%s_to_%s.out" % gpair)
    logerr = os.path.join(outdir, "%s_to_%s.err" % gpair)
    with open(logout, "w") as fo, open(logerr, "w") as fe:
        return (thiscmd, subprocess.Popen(thiscmd, stdout = fo, stderr = fe))

procs = map(make_proc, compare_groups)

anyerr = False

for (cmd, p) in procs:
    if p.wait() != 0:
        print >>sys.stderr, "Warning: proc %s failed" % cmd

