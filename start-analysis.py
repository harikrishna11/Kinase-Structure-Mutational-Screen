#!/usr/bin/env python

import sys
import os.path
import subprocess
import checkjobfile
import _mysql
import tempfile
import os

## This is the master script for starting an NSMD analysis run. This stage runs after the main GROMACS
## run has completed, and requires also that a corresponding wild-type run has completed (if this run *is*
## wild-type then the analysis process is self-sufficient)
## It takes care of checking that the relevant runs have indeed completed, checking that the kinase is sufficiently
## structurally annotated in the database for analysis, and allowing the user to supply structural annotations if missing.

nsmd_root = "/data/snc/nstephenson/nsmd"
nsmd_makefile = os.path.join(nsmd_root, "scripts/Makefile")

# The database is needed to check for structural annotations (e.g. that the DFG motif has been located).
db = _mysql.connect(host = "acbbdb1.picr.man.ac.uk", user = "nsmd", passwd = "nsmdP123", db = "nsmd")

# This argument structure is the same as start-nsmd.py. See the comments on that script for the PDB:CHAIN and mutation syntax also.
if len(sys.argv) < 4:
    print >>sys.stderr, "Usage: start-analysis.py gene active|inactive {mutation1,mutation2... | wild} [PDB:CHAIN]"
    sys.exit(1)

if sys.argv[2] == "active":
    use_active_model = True
elif sys.argv[2] == "inactive":
    use_active_model = False
else:
    raise Exception("Argument 2 must be 'active' or 'inactive'")

if len(sys.argv) > 4:
    pdbforced = True
    pdbid, pdbchain = sys.argv[4].split(":")
else:
    pdbforced = False
    pdbid, pdbchain = None, None

runsuffix = "" if use_active_model else "_inactive"
if pdbforced:
    runsuffix = runsuffix + "_%s_%s" % (pdbid, pdbchain)

gene = sys.argv[1]
mutation = sys.argv[3]

gene_bits = gene.split("-")
if len(gene_bits) == 1:
    where_clause = "gene = '%s'" % gene_bits[0]
elif len(gene_bits) == 2:
    where_clause = "gene = '%s' and offset_canonical = '%s'" % tuple(gene_bits)
else:
    raise Exception("Bad syntax: gene can contain at most one - character")

# Fetch any motif annotation to determine whether we've annotated this kinase at all yet:

query = 'select dfg_motif_offset,dxxxxg_motif_offset from genes join kinases on genes.id = gene_id where ' + where_clause
db.query(query)

rows = db.store_result().fetch_row(maxrows=0)

if len(rows) == 0:
    raise Exception("No such gene %s" % gene)
elif len(rows) > 1:
    print >>sys.stderr, "Specified gene has more than one kinase domain. You should specify which offset to use by replacing %s with one of the following:" % gene
    for r in rows:
        print >>sys.stderr, "%s-%s" % (gene, r[1])
    sys.exit(1)

if rows[0][0] is None:
    print >>sys.stderr, "We don't have feature annotation for this kinase in the database yet. Trying to automatically add it..."

if rows[0][1] is None:
    print >>sys.stderr, "dxxxxg_motif_offset entry is missing in the database. Trying to automatically add it..." 

if None in rows[0]:

    # This script tries to automatically find all the interesting motifs (DFG, HRD, ...) by looking for unique best matches.
    # However often it will (partially) fail and appeal to the user for help resolving the ambiguity.

    script = os.path.join(os.path.dirname(sys.argv[0]), "report-all-kinase-features.py")
    active_arg = "manual" if pdbforced else sys.argv[2]
    genespec = ("%s:%s:%s" % (gene_bits[0], pdbid, pdbchain)) if pdbforced else gene
    tempdir = tempfile.mkdtemp()
    outcsv = os.path.join(tempdir, "report.csv")
    errlog = os.path.join(tempdir, "stderr.log")
    cmdline = [script, active_arg, "/tmp", genespec]
    with open(outcsv, "w") as of, open(errlog, "w") as ef:
        subprocess.check_call(cmdline, stdout = of, stderr = ef)
    with open(outcsv, "r") as f:
        content = f.read()
        nlines = len([l for l in content.split("\n") if len(l.strip()) > 0])
        if nlines <= 1:
            print >>sys.stderr, "Search command", " ".join(["'%s'" % s for s in cmdline]), "did not return any records"
            sys.exit(1)
        elif content.find("Fail") != -1:
            print >>sys.stderr, "Couldn't automatically find all needed features. Will print partial results; please manually identify the missing features:"
            print >>sys.stderr, content
            sys.exit(1)
        else:
            print >>sys.stderr, "Success! Inserting the following records into the database:"
            print >>sys.stderr, content

    # Write the new kinase annotation to the database:

    script = os.path.join(os.path.dirname(sys.argv[0]), "import-kinase-feature-spreadsheet.py")
    cmdline = [script, outcsv, "auto"]
    subprocess.check_call(cmdline)

# We need the corresponding wild-type directory (which might be equal to result_dir if this *is* a wild-type run) for some mutant-vs-wild analysis.

result_dir = os.path.join(nsmd_root, "results/%s/%s%s" % (gene, mutation, runsuffix))
wild_dir = os.path.join(nsmd_root, "results/%s/wild%s" % (gene, runsuffix))

# These two files are hallmarks of a finished analysis run. NPT_CA.gro is separately checked-for because it was introduced
# after some runs were already complete, so some runs may need re-execution to produce it despite having notionally completed.

need_files = ["main_sim_400ns_done.flag", "NPT_CA.gro"]

# Check that both the target run and the corresponding wild run are ready for analysis.

for d in [result_dir, wild_dir]:
    for f in need_files:
        needed = os.path.join(d, f)
        if not os.path.exists(needed):
            print >>sys.stderr, "%s not found (indicates needed analysis not yet finished)" % needed
            sys.exit(1)

# The wild-type structure needs also to have been *analysed*, not just run.

if mutation != "wild":
    needed = os.path.join(wild_dir, "400ns_specific_graphs.done")
    if not os.path.exists(needed):
        print >>sys.stderr, "%s does not exist. Run wild type analysis before any mutants." % needed
        sys.exit(1)

# Start a cluster job to do the actual work:

jobname = "nsmd_preemptable.analysis.%s.%s%s" % (gene, mutation, runsuffix)
outpath = os.path.join(result_dir, "analysis.out")
errpath = os.path.join(result_dir, "analysis.err")

command = ["/bin/bash", "-c", "echo make -f %s -C %s 400ns_analyses WILDDIR=%s | qsub -h -N %s -o %s -e %s -l nodes=1:ppn=16 -m a -M Natalie.Stephenson@cruk.manchester.ac.uk -lwalltime=20:00:00:00" % (nsmd_makefile, result_dir, wild_dir, jobname, outpath, errpath)]

if os.getenv("DRYRUN") is not None:
    print command
    sys.exit(0)

jobfile = os.path.join(result_dir, "analysis_jobid.txt")
checkjobfile.checkjobfile(jobfile)

try:
    with open(jobfile, "w") as f:
        subprocess.check_call(command, stdout = f)
except Exception as e:
    try:
        os.unlink(jobfile)
    except:
        pass
    raise e



