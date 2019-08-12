#!/apps/scripts/nsmd/deps/modeller-9.14-install/bin/modpy.sh /apps/modules/pkg/apps/python/2.7.6/gcc-4.4.7/bin/python

import sys
import modeller

if len(sys.argv) < 3 or len(sys.argv) > 5:
    print >>sys.stderr, "Usage: compare-pdbs.py 1.pdb 2.pdb [1chain [2chain]]"
    sys.exit(1)

while len(sys.argv) < 5:
    sys.argv.append(None)

env = modeller.environ()
aln = modeller.alignment(env)

def file_to_model(pdbfile, chain):
    if chain is None:
        return modeller.model(env, file=pdbfile)
    else:
        return modeller.model(env, file=pdbfile, model_segment=('FIRST:%s' % chain, 'LAST:%s' % chain))

mdl1 = file_to_model(sys.argv[1], sys.argv[3])
mdl2 = file_to_model(sys.argv[2], sys.argv[4])

aln.append_model(mdl1, align_codes='aln1', atom_files=sys.argv[1])
aln.append_model(mdl2, align_codes='aln2', atom_files=sys.argv[2])
aln.align()
aln.compare_structures(compare_mode = 1)

