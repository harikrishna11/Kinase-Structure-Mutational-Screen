#!/usr/bin/python

import sys
import _mysql
import pdbquery
import mutate
import tempfile
import shutil

## This script is used to populate the various PDB-related columns in the 'kinases' and 'inactive_pdbs' tables.
## It searches for PDB models relating to the right gene using PDB's web API, then downloads the PDBs recommended,
## checks that they overlap with the right part of the protein sequence (the kinase domain), and categorises them
## as active, inactive or uncertain configurations.

## The active/inactive classification uses the DFG motif, which is currently internally discovered because this step usually
## comes before the motif-annotation step. This could potentially be fixed if it becomes sufficiently annoying.

## In practice it is currently common for Natalie to specify manual PDB entries instead of relying on this script's findings,
## so it may be less relevant than in the past.

db = _mysql.connect(host = "acbbdb1.picr.man.ac.uk", user = "nsmd", passwd = "nsmdP123", db = "nsmd")

if sys.argv[1] == "active":
        find_active_model = True
elif sys.argv[1] == "inactive":
        find_active_model = False
else:
        raise Exception("Must specify 'active' or 'inactive'")

# Specifying a work directory is recommended for manual debugging, since the second run will re-used already-fetched PDB files.

if len(sys.argv) >= 3:
        workdirtemp = False
        workdir = sys.argv[2]
else:
        workdirtemp = True
        workdir = tempfile.mkdtemp()

def strornull(f):
        if f is None:
                return "NULL"
        else:
                return str(f)

update_table = "kinases" if find_active_model else "inactive_pdbs"

# The bulk of the actual work is done by pdbquery.getpdbinfo and pdbquery.getbestentry. See that file for further commentary.

def create_gene(recordid, name, pdbdict, seq, domain, uniprotid):

        if len(pdbdict) == 0:
                print >>sys.stderr, "No PDB entries for", name
                bestpdb = None
        else:
                bestpdb = pdbquery.getbestentry(name, pdbdict, seq, domain, workdir, uniprotid, find_active_model = find_active_model)

        if bestpdb is None:
                db.query("update %s set pdbid = 'NOTFOUND' where id = %d" % (update_table, recordid))
        else:
                db.query("update %s set pdbid = '%s', pdbchain = '%s', pdbdesc = '%s', pdbres = %g, pdbkinasecoverage = %g, pdbrvalue = %s, pdbrfree = %s where id = %d" % (update_table, bestpdb["id"], bestpdb["best_chain"], db.escape_string(bestpdb["description"]), bestpdb["resolution"], bestpdb["coverage"], strornull(bestpdb["rvalue"]), strornull(bestpdb["rfree"]), recordid))
        db.store_result()

# Try to find a PDB record for all genes/kinase domains in the database that don't have one already.

db.query("select kinases.id, gene, uniprotid, canonical_sequence, offset_canonical, length, %s.pdbid from genes inner join kinases on genes.id = kinases.gene_id inner join inactive_pdbs on kinases.id = inactive_pdbs.id where offset_canonical is not null" % update_table)
for recordid, gene, uniprotid, seq, kinase_offset, kinase_len, pdbid in db.store_result().fetch_row(maxrows=0):

        if pdbid is not None:
                continue

        kinase_offset = int(kinase_offset)
        recordid = int(recordid)
        kinase_len = int(kinase_len)

        domain = seq[kinase_offset : kinase_offset + kinase_len]

        pdbresults = pdbquery.getpdbinfo([uniprotid])

        # Result is a dict of dicts: uniprotid -> pdbid -> result descriptions
        pdbdict = pdbresults[uniprotid]
        create_gene(recordid, gene, pdbdict, seq, domain, uniprotid)

if workdirtemp:
        shutil.rmtree(workdir)
