#!/usr/bin/env python

import sys
import os.path
import shutil

if len(sys.argv) < 3:
    print >>sys.stderr, "Usage: nsmd-export.py file1 [file2 file3 ...] destination"
    sys.exit(1)

nsmd_result_root = "/mnt/gpfs2/data/snc/nstephenson/nsmd/results/"

target = sys.argv[-1]
sources = sys.argv[1:-1]

if not os.path.isdir(target):
    print >>sys.stderr, "%s: not a directory" % target
    sys.exit(1)

for f in sources:
    if not os.path.exists(f):
        print >>sys.stderr, "%s: not found" % f
        sys.exit(1)

for f in sources:
    abspath = os.path.realpath(os.path.abspath(f))
    if not abspath.startswith(nsmd_result_root):
        print >>sys.stderr, "Warning: %s not in nsmd results directory; copying with existing name" % f
        targetname = os.path.basename(f)
    else:
        relpath = abspath[len(nsmd_result_root):]
        bits = relpath.split("/")
        targetname = "_".join(bits)
    shutil.copyfile(f, os.path.join(target, targetname))
        
