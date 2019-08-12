#!/usr/bin/env python

import sys
import _mysql
import subprocess
import tempfile
import os.path
import modelutils

if len(sys.argv) >= 2:
    workdir = sys.argv[1]
    sys.argv = sys.argv[1:]
    workdirtemp = False
else:
    workdir = tempfile.mkdtemp()
    workdirtemp = True

mydir = os.path.dirname(os.path.realpath(__file__))

commit = False

if len(sys.argv) == 2:
    if sys.argv[1] == "commit":
        commit = True
        print >>sys.stderr, "Will commit changes"
    else:
        raise Exception("Usage: check-bridge-distances.py [workdir] [commit]")

db = _mysql.connect(host = "acbbdb1.picr.man.ac.uk", user = "nsmd", passwd = "nsmdP123", db = "nsmd")

query = 'select kinases.id, gene, pdbid, pdbchain, offset_canonical, bridge_glut_offset, bridge_lys_offset, bridge_distance from genes join kinases on genes.id = gene_id where pdbid is not null and pdbid != "NOTFOUND" and bridge_glut_offset is not null and bridge_lys_offset is not null'

devnull = open("/dev/null", "w")

db.query(query)

for recordid, gene, pdbid, pdbchain, offset, bridge_glut, bridge_lys, bridge_distance in db.store_result().fetch_row(maxrows=0):

    try:
        if bridge_distance is not None:
            bridge_distance = float(bridge_distance)
    except Exception as e:
        print >>sys.stderr, "Can't convert", bridge_distance
        raise e

    pdbfile = modelutils.get_pdb_file(pdbid, workdir)

    params = [os.path.join(mydir, "base-distance.py"), pdbfile, pdbchain, bridge_glut, bridge_lys]
    proc = subprocess.Popen(params, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (outs, errs) = proc.communicate()
    
    try:
        new_distance = float(outs)

        if bridge_distance is None or abs(new_distance - bridge_distance) > 0.001:
            print "Deviation: %s/%s/%s/%s %s-%s had %s but now measures %g" % (gene, offset, pdbid, pdbchain, bridge_glut, bridge_lys, bridge_distance, new_distance)
            if commit:
                db.query("update kinases set bridge_distance = %g where id = %d" % (new_distance, int(recordid)))
                db.store_result()

    except ValueError:

        print >>sys.stderr, "Call", params, "failed"
        print >>sys.stderr, "Out:", outs
        print >>sys.stderr, "Err:", errs

if workdirtemp:
    shutil.rmtree(workdir)
