#!/usr/bin/env python

import _mysql
import subprocess
import os
import os.path
import sys

rootdir = "/data/snc/nstephenson/nsmd/results"

db = _mysql.connect(host = "acbbdb1", user = "nsmd", passwd = "nsmdP123", db = "nsmd")

mydir = os.path.realpath(__file__)

for gene in os.listdir(rootdir):

    if gene == "ABL1_MODELS":
        continue
    
    genebits = gene.split("-")
    if len(genebits) < 1 or len(genebits) > 2:
        raise Exception("Unexpected gene syntax %s" % gene)

    where_clause = "gene = '%s'" % genebits[0]
    if len(genebits) == 2:
        where_clause = "%s and offset_canonical = %s" % (where_clause, genebits[1])

    db.query("select offset_canonical from genes join kinases on gene_id = genes.id where %s" % where_clause)
    try:
        offset = int(db.store_result().fetch_row(maxrows=0)[0][0])
    except IndexError:
        print >>sys.stderr, "No gene", gene

    for runspec in os.listdir(os.path.join(rootdir, gene)):

        runpath = os.path.join(rootdir, gene, runspec)
        jobfile = os.path.join(runpath, "jobid.txt")
        donefile = os.path.join(runpath, "main_sim_400ns_done.flag")
        if (not os.path.exists(jobfile)) or os.path.exists(donefile):
            continue

        with open(jobfile, "r") as f:
            oldjobid = f.read().strip()
            with open("/dev/null", "w") as nullout:
                checkret = subprocess.Popen(["/apps/modules/pkg/clusterware/torque/5.1.0-1/gcc-4.4.7/bin/qstat", oldjobid], stdout=nullout, stderr=subprocess.STDOUT).wait()
            if checkret == 0:
                print >>sys.stderr, "Skipped %s/%s already in progress" % (gene, runspec)
                continue

        restart_params = ["/bin/bash", os.path.join(mydir, "start-nsmd.sh"), gene]
        runspec_bits = runspec.split("_")
        mutspec = runspec_bits[0]
        runspec_bits = runspec_bits[1:]
        
        if len(runspec_bits) > 0 and runspec_bits[0] == "inactive":
            restart_params.append("inactive")
            runspec_bits = runspec_bits[1:]
        else:
            restart_params.append("active")
            
        def add_offset(mut, offset):
            if mut == "wild":
                return mut
            oldoff = int(mut[1:-1])
            return "%s%d%s" % (mut[0], oldoff + offset, mut[-1])

        mutspec = ",".join([add_offset(m, offset) for m in mutspec.split(",")])

        restart_params.append(mutspec)
        
        if len(runspec_bits) > 0:
            restart_params.append("%s:%s" % tuple(runspec_bits))

        print >>sys.stderr, "Issuing restart:", restart_params
        try:
            subprocess.check_call(restart_params)
        except Exception as e:
            print e


            
        
