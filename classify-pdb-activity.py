#!/usr/bin/env python

import pdbquery
import shutil
import tempfile
import os.path
import sys
import subprocess

if len(sys.argv) < 4:
    print >>sys.stderr, "Usage: classify-pdb-activity.py uniprot_code kinase_domain_start kinase_domain_end"
    sys.exit(1)

dom_start = int(sys.argv[2])
dom_end = int(sys.argv[3])

workdir = tempfile.mkdtemp()

results = dict()

def overlap(s1, e1, s2, e2):
    if e1 < s2:
        return False
    elif s1 > e2:
        return False
    return True

pdbs = pdbquery.getpdbinfo([sys.argv[1]], verbose=True)[sys.argv[1]]
for (pdb, alignments) in pdbs.iteritems():

    for details in alignments:

        chain_coverage = []

        print details["chains"]

        for (chain, regions) in details["chains"].iteritems():
            for (start, end) in regions:
                if overlap(dom_start, dom_end, start, end):
                    if len(chain_coverage) == 0 or chain_coverage[-1] != chain:
                        chain_coverage.append(chain)

        if len(chain_coverage) == 0:
            print >>sys.stderr, "Ignore", pdb, "with no kinase domain coverage"
            continue
        elif len(chain_coverage) != 1:
            print >>sys.stderr, "Warning!", pdb, "has multiple chains mentioned within one alignment"
        else:
            print >>sys.stderr, "Checking", pdb

        pdbdata = pdbquery.get("www.rcsb.org", "/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=%s" % pdb)
        pdbfname = os.path.join(workdir, "%s.pdb" % pdb)
        with open(pdbfname, "w") as f:
            f.write(pdbdata)

        for chain in chain_coverage:

            p = subprocess.Popen(["./check-active-dfg.py", pdbfname, chain], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (outs, errs) = p.communicate()
            distance_report = ", ".join([l for l in errs.split("\n") if l.find("Distance") != -1])
            ret = p.returncode
            if ret == 0:
                result = "Active"
            elif ret == 1:
                result = "Inactive"
            else:
                result = "Error"
            if result not in results:
                results[result] = []
            results[result].append((pdb, chain, distance_report))

for (result, pdbchains) in results.iteritems():
    print result
    for (pdb, chain, distance_report) in pdbchains:
        print "\t%s.%s (%s)" % (pdb, chain, distance_report)

shutil.rmtree(workdir)
