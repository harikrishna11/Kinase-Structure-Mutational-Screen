#!/usr/bin/env python

import subprocess
import sys
import _mysql
import pexpect
import read_groups
import parse_genespec

if len(sys.argv) <= 2:
    print >>sys.stderr, "Usage: make_gro_index.py gene rmsd|dist make_ndx -f sim.gro -o index.ndx ..."
    sys.exit(1)

where_clause = parse_genespec.get_whereclause(sys.argv[1])

db = _mysql.connect(host = "acbbdb1", user = "nsmd", passwd = "nsmdP123", db = "nsmd")
db.query("select offset_canonical, ape_motif_offset, dfg_motif_offset, hrd_motif_offset, gxgxxg_motif_offset, vaik_motif_offset, activation_loop_start, activation_loop_end, achelix_start, achelix_end, bridge_glut_offset, bridge_lys_offset, dxxxxg_motif_offset, rspine_start_offset, rspine_ac_offset from genes join kinases on genes.id = gene_id where " + where_clause)

rows = db.store_result().fetch_row(maxrows=0)

if len(rows) == 0:
    raise Exception("No such gene %s" % sys.argv[1])
elif len(rows) > 1:
    raise Exception("Gene %s is ambiguous; specify a kinase domain offset" % sys.argv[1])

kinase_offset = int(rows[0][0])

ape_off, dfg_off, hrd_off, gxgxxg_off, vaik_off, actloop_start, actloop_end, achelix_start, achelix_end, bridge_glut, bridge_lys, dxxxxg_off, rs1_off, rs2_off = map(lambda x : int(x) - kinase_offset, rows[0][1:])

if sys.argv[2] == "rmsd":
    groups = read_groups.get_rmsd_groups(ape_off, dfg_off, hrd_off, gxgxxg_off, vaik_off, actloop_start, actloop_end, achelix_start, achelix_end, bridge_glut, bridge_lys)
elif sys.argv[2] == "all":
    groups = read_groups.get_dist_groups(ape_off, dfg_off, hrd_off, gxgxxg_off, vaik_off, actloop_start, actloop_end, achelix_start, achelix_end, bridge_glut, bridge_lys, dxxxxg_off, rs1_off, rs2_off)
    groups = groups + read_groups.get_rmsd_groups(ape_off, dfg_off, hrd_off, gxgxxg_off, vaik_off, actloop_start, actloop_end, achelix_start, achelix_end, bridge_glut, bridge_lys)
else:
    raise Exception("Argument 2 must be rmsd or all")

idx_proc = pexpect.spawn(" ".join(sys.argv[3:]))

idx_proc.expect("\r\n>")

for (nm, defn) in groups:

    print "Creating", nm, defn

    idx_proc.sendline("%s" % defn)
    idx_proc.expect("\r\n>")

    groupnum = None
    for l in idx_proc.before.split("\n"):
        if l.find("atoms") != -1:
            bits = l.split()
            groupnum = bits[0]

    if groupnum is not None:
        idx_proc.sendline("name %s %s" % (groupnum, nm))
    else:
        print >>sys.stderr, "Couldn't name", nm
    idx_proc.expect("\r\n>")

idx_proc.sendline("q")
idx_proc.wait()
    
    


                

