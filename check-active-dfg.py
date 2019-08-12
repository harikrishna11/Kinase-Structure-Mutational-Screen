#!/apps/scripts/nsmd/deps/modeller-9.14-install/bin/modpy.sh python

import modeller
import sys
import math
import modelutils

## This checks whether a given PDB file / chain apears to be an active / inactive kinase model.
## It uses a simple naive search for the DFG and HRD motifs -- this should probably be replaced with
## reference to the database, as it might even find motifs that are outside the kinase domain!

# Exit codes:
# 0 -> active conformation
# 1 -> inactive
# other -> error

if len(sys.argv) < 3:
    print >>sys.stderr, "Usage: check-active-dfg.py in.pdb chain [distance]"
    sys.exit(2)

pdbcode = sys.argv[1]
chain = sys.argv[2]

print >>sys.stderr, "Checking %s:%s" % (pdbcode, chain)

env = modeller.environ()
mdl = modeller.model(env, file = pdbcode, model_segment = ('FIRST:' + chain, 'LAST:' + chain))

seq = "".join([res.code for res in mdl.residues])

def find_or_die(haystack, needle):
    off = haystack.find(needle)
    if off == -1:
        print >>sys.stderr, "Sequence", needle, "not found"
        sys.exit(2)
    else:
        print >>sys.stderr, "Found", needle, "at offset", off
    return off

dfg_off = find_or_die(seq, "DFG")
hrd_off = find_or_die(seq, "HRD")

tests = [
#    (dfg_off, "OD1", dfg_off + 2, "N"),
    (dfg_off + 1, "O", dfg_off + 4, "N"),
    (hrd_off + 1, "NE", dfg_off + 3, "O"),
#    (hrd_off + 1, "NH2", dfg_off + 3, "O")
]

dists = [modelutils.atom_dist(mdl, *t) for t in tests]
if len(sys.argv) >= 4 and sys.argv[3] == "distance":
    print >>sys.stderr, max(dists)
elif max(dists) < 4:
    print >>sys.stderr, "Active"
    sys.exit(0)
elif max(dists) > 6:
    print >>sys.stderr, "Inactive"
    sys.exit(0)
else:
    print >>sys.stderr, "Unknown"
    sys.exit(2)
