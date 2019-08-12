#!/apps/scripts/nsmd/deps/modeller-9.14-install/bin/modpy.sh python

import modeller
import modelutils
import sys

# Measure the distance between two bases in a given PDB file/chain. The given offsets are
# 1-based against the Uniprot structure and are adjusted to PDB coordinates first.

if len(sys.argv) != 5:
    print >>sys.stderr, "Usage: base-distance.py structure.pdb chain offset1 offset2"
    sys.exit(1)

pdbfile = sys.argv[1]
chain = sys.argv[2]
off1 = int(sys.argv[3])
off2 = int(sys.argv[4])

pdb_to_uniprot = modelutils.read_pdb_to_uniprot(pdbfile, chain)
uniprot_to_pdb = dict([(v, k) for (k, v) in pdb_to_uniprot.iteritems()])

off1 = uniprot_to_pdb[off1]
off2 = uniprot_to_pdb[off2]

# Supress verbose version notice
with open("/dev/null", "w") as fnull:
    oldout = sys.stdout
    sys.stdout = fnull
    env = modeller.environ()
    mdl = modeller.model(env, file = pdbfile, model_segment = ('FIRST:' + chain, 'LAST:' + chain))
    sys.stdout = oldout

def find_res(off):

    match = [res for res in mdl.residues if int(res.num) == off]
    if len(match) != 1:
        raise Exception("Found %d residues with PDB offset %d" % (len(match), off))

    return match[0]

res1 = find_res(off1)
res2 = find_res(off2)

print modelutils.closest_res_dist(res1, res2)

