#!/usr/bin/env python

import sys
import os
import os.path
import subprocess
import _mysql
import mutate
import checkjobfile
import parse_genespec

## This is the master script for kicking off a new NSMD run. It takes care of ensuring duplicate runs aren't created, prepares a run directory under /data/snc/nstephenson/results,
## and finally queues a cluster job that will run 'make' with the run directory as the cwd.

nsmd_root = "/data/snc/nstephenson/nsmd"
nsmd_makefile = os.path.join(nsmd_root, "scripts/Makefile")

if len(sys.argv) < 4 or len(sys.argv) > 5:
    print >>sys.stderr, "Usage: start-nsmd.py gene active|inactive {mutation1,mutation2... | wild} [PDB:CHAIN]"
    sys.exit(1)

# The DB is used to determine where the kinase domain sits within the whole protein sequence, and to find the default PDB model for this kinase if the user hasn't given one manually.
db = _mysql.connect(host = "acbbdb1", user = "nsmd", passwd = "nsmdP123", db = "nsmd")

gene = sys.argv[1]
if sys.argv[2] == "active":
    use_active_model = True
elif sys.argv[2] == "inactive":
    use_active_model = False
else:
    raise Exception("Argument #2 must be 'active' or 'inactive'")

mutation = sys.argv[3]

# 'gene' might just be a gene name, or might have the form ABC1-NNN, where NNN is an integer offset into the protein.
# The latter form is needed for proteins with multiple kinase domains, like the Janus kinases. get_whereclause returns an appropriate SQL query to select the row in question.
where_clause = parse_genespec.get_whereclause(gene)

# Preferred PDB models for an active kinase configuration are stored directly in the 'kinases' table.
pdbtable = "kinases" if use_active_model else "inactive_pdbs";

db.query("select canonical_sequence, offset_canonical, length, %s.pdbid, %s.pdbchain from genes join kinases on gene_id = genes.id join inactive_pdbs on kinases.id = inactive_pdbs.id where %s" % (pdbtable, pdbtable, where_clause))
rows = db.store_result().fetch_row(maxrows=0)

if len(rows) == 0:
    raise Exception("No such gene %s" % gene)
elif len(rows) > 1:
    print >>sys.stderr, "Specified gene has more than one kinase domain. You should specify which offset to use by replacing %s with one of the following:" % gene
    for r in rows:
        print >>sys.stderr, "%s-%s" % (gene, r[1])
    sys.exit(1)

sequence, offset, length, pdbid, pdbchain = rows[0]
pdbforced = False

# Did the user specify a model to use instead of the preferred one stored in the DB?
if len(sys.argv) == 5:
    pdbbits = sys.argv[4].split(":")
    if len(pdbbits) != 2:
        print >>sys.stderr, "Forced PDB specification must have form ID:CHAIN"
        sys.exit(1)
    pdbid, pdbchain = pdbbits
    pdbforced = True

if pdbid == "NOTFOUND" or pdbid == "" or pdbid is None:
    raise Exception("We don't seem to have a candidate PDB model for gene %s" % gene)

# Subsequence kinase domain:
sequence = sequence[int(offset) : int(offset) + int(length)]

muts = mutation.split(",")

# The special mutation-spec 'wild' means no mutations (wild type).
if not (len(muts) == 1 and muts[0] == "wild"):

    # Rewrite mutations relative to the kinase domain. Mutations given on the command-line are specified protein-relative; here we'll rewrite them
    # to be kinase-domain-relative, and check they actually agree with the given sequence (mutations are specified like ANNNB, where A and B
    # are the old and new residues respectively, and NNN is an offset.

    parsed_muts = map(mutate.parsespec, muts)

    for oldbase, newbase, mutoff in parsed_muts:
        if mutoff - 1 < int(offset) or mutoff - 1 >= int(offset) + int(length):
            print >>sys.stderr, "Given mutation %s%d%s out of range for the kinase domain, which spans %d-%d" % (oldbase, mutoff, newbase, int(offset) + 1, int(offset) + int(length) + 1)
            sys.exit(1)

    offset_muts = [(oldbase, newbase, mutoff - int(offset)) for (oldbase, newbase, mutoff) in parsed_muts]

    for oldbase, newbase, mutoff in offset_muts:
        if sequence[mutoff - 1] != oldbase:
            raise Exception("Given mutation %s%d%s disagrees with base sequence, which has %s at that position" % (oldbase, mutoff + int(offset), newbase, sequence[mutoff - 1]))

    mutation = ",".join(["%s%d%s" % (oldbase, mutoff, newbase) for (oldbase, newbase, mutoff) in offset_muts])

# Create a run directory, named after the mutation spec, followed by _inactive if investigating an inactive kinase, and/or _PDBID:PDBCHAIN if
# a manual PDB ID has been supplied.

runsuffix = "" if use_active_model else "_inactive"
if pdbforced:
    runsuffix = runsuffix + "_%s_%s" % (pdbid, pdbchain)

runpath = os.path.join(nsmd_root, "results/%s/%s%s" % (gene, mutation, runsuffix))
needed_files = ["base_sequence.txt", "base_sequence_mutation.txt", "base_model.pdb", "base_model_chain.txt"]
needed_files_abs = [os.path.join(runpath, x) for x in needed_files]

# This is for tracking the run's cluster job, for resuming after a run dies for some reason.
jobfile = os.path.join(runpath, "jobid.txt")

if os.path.exists(runpath):

    print >>sys.stderr, "Warning: result directory %s already exists" % runpath
    if not all(map(os.path.exists, needed_files_abs)):
        print >>sys.stderr, "At least one of %s does not exist. Most likely the run directory was part-created; recommend deleting it and starting over." % needed_files
        sys.exit(1)

    checkjobfile.checkjobfile(jobfile)

else:

    # Write the basic information needed by the Makefile to start a new run.

    os.makedirs(runpath)

    seqfile, mutfile, pdbfile, chainfile = needed_files_abs

    with open(seqfile, "w") as f:
        f.write(sequence)

    with open(mutfile, "w") as f:
        if mutation != "wild":
            f.write(mutation)

    print >>sys.stderr, "Fetch PDB file %s" % pdbid
    subprocess.check_call(["/usr/bin/wget", "-O", pdbfile, "http://www.rcsb.org/pdb/files/%s.pdb" % pdbid])

    with open(chainfile, "w") as f:
        f.write(pdbchain)

# Start the cluster job.

jobname = "nsmd_preemptable.%s.%s%s" % (gene, mutation, runsuffix)

outpath = os.path.join(runpath, "nsmd.stdout.log")
errpath = os.path.join(runpath, "nsmd.stdout.log")

command = ["/bin/bash", "-c", "echo make -f %s -C %s | qsub -h -N %s -o %s -e %s -l nodes=1:ppn=16 -m a -M Natalie.Stephenson@cruk.manchester.ac.uk -lwalltime=20:00:00:00" % (nsmd_makefile, runpath, jobname, outpath, errpath)]

try:
    with open(jobfile, "w") as f:
        subprocess.check_call(command, stdout = f)
except Exception as e:
    try:
        os.unlink(jobfile)
    except:
        pass
    raise e
