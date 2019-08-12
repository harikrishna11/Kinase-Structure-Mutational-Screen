#!/apps/scripts/nsmd/deps/modeller-9.14-install/bin/modpy.sh /apps/modules/pkg/apps/python/2.7.6/gcc-4.4.7/bin/python

import sys
import tempfile
import os
import os.path
import mutate
import modeller
import modeller.automodel
import subprocess
import shutil

if len(sys.argv) < 5 or len(sys.argv) > 6 or not sys.argv[2].endswith(".pdb"):
    print >>sys.stderr, "Usage: create-model.py in_sequence model.pdb model_chain out.pdb [mutation]"
    sys.exit(1)

mydir = os.path.dirname(os.path.realpath(__file__))

# The os.chdir later will cause trouble if the .pdb file path is not absolute.
sys.argv[2] = os.path.abspath(sys.argv[2])

td = tempfile.mkdtemp()

if len(sys.argv) == 6 and len(sys.argv[5]) != 0:
    for mut in sys.argv[5].split(","):
        old_bases, new_bases, base_pos = mutate.parsespec(mut)
        sys.argv[1] = mutate.applymutation(sys.argv[1], old_bases, new_bases, base_pos)
    
pir_file = os.path.join(td, "seq.pir")
with open(pir_file, "w") as f:
    f.write(">P1;INPUT\n")
    f.write("sequence:INPUT:::::::0.00: 0.00\n")
    f.write(sys.argv[1])
    f.write("*\n")

align_code = os.path.split(sys.argv[2])[1][:-4] + sys.argv[3]

env = modeller.environ()
aln = modeller.alignment(env)
mdl = modeller.model(env, file=sys.argv[2], model_segment=('FIRST:%s' % sys.argv[3], 'LAST:%s' % sys.argv[3]))

align_file = os.path.join(td, "seqmodel.ali")

aln.append_model(mdl, align_codes=align_code, atom_files=sys.argv[2])
aln.append(file=pir_file, align_codes='INPUT')
aln.align2d()
aln.write(file=align_file, alignment_format='PIR')

# Automodel spits crap out into the cwd :\
script_dir = os.getcwd()
compare_script = os.path.join(mydir, "compare-pdbs.py")
out_pdb_path = os.path.join(script_dir, sys.argv[4])
os.chdir(td)

a = modeller.automodel.automodel(env, alnfile=align_file, knowns=align_code, sequence='INPUT')
a.starting_model = 1
a.ending_model = 10
a.make()

def output_model_rmsd(out_mdl):

    compare_out = subprocess.check_output([compare_script, out_mdl["name"], sys.argv[2]])
    lines = [l.strip() for l in compare_out.split("\n")]
    found_header = False
    for l in lines:
        if l.find("Position comparison (FIT_ATOMS)") != -1:
            found_header = True
        elif found_header and l.startswith("aln1") and len(l.split()) == 3:
            bits = l.split()
            return float(bits[2])

models_with_rmsd = [(m, output_model_rmsd(m)) for m in a.outputs]
models_with_rmsd = sorted(models_with_rmsd, key = lambda x : x[1])
best_mod = models_with_rmsd[0]

print >>sys.stderr, "Using model", best_mod[0]["name"], "with RMSD", best_mod[1]
shutil.copyfile(best_mod[0]["name"], out_pdb_path)
    
shutil.rmtree(td)
