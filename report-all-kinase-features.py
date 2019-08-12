#!/usr/bin/env python

import _mysql
import os
import subprocess
import json
import shutil
import csv
import sys
import tempfile
import os.path

## Given a kinase specification, try to find interesting features (e.g. the DFG motif). This is used from start-analysis.py if a protein is not already annotated.
## It can also usefully be run manually in order to pre-populate features after new kinases have been imported into the database.

# Running with argument 'active' or 'inactive' will 

mydir = os.path.realpath(__file__)

if sys.argv[1] == "active":
    pdbtable = "kinases"
    have_pdbtable = True
elif sys.argv[1] == "inactive":
    pdbtable = "inactive_pdbs"
    have_pdbtable = True
elif sys.argv[1] == "manual":
    pdbtable = "kinases" # Harmless default, not used
    have_pdbtable = False
else:
    raise Exception("First argument must be 'active', 'inactive' or 'manual'")

workdir = sys.argv[2]
workdirtemp = False

# If either 'active' or 'inactive' are specified then we'll try to find features for all kinases that have a PDB model in the database.
# Otherwise the user must manually specify arguments of the form gene-name:pdbid:pdbchain.
if len(sys.argv) == 3 and not have_pdbtable:
    raise Exception("In manual mode gene specs must be given as arguments")

db = _mysql.connect(host = "acbbdb1.picr.man.ac.uk", user = "nsmd", passwd = "nsmdP123", db = "nsmd")

workfile = os.path.join(workdir, "out.json")

# CSV columns that we will output. Intentionally user-friendly rather than matching the database column names so this can be directly handed to the user for completion if any entries are missing.
cols = ["Gene", "Kinase offset", "PDB ID", "APE Start", "DFG Start", "HRD Start", "GxGxxG Start", "VAIK Start", "VAIK Sequence", "Activation Loop Start", "Activation Loop End", "aC Helix Start", "aC Helix End", "Bridge Glutamine", "Bridge Lysine", "Bridge Closest Atoms", "DxxxxG Start", "R-Spine B4", "R-Spine aC"]

devnull = open("/dev/null", "w")
writer = csv.DictWriter(sys.stdout, fieldnames=cols)
writer.writeheader()

def parse_genespec(spec):
    bits = spec.split(":")
    if len(bits) != 1 and len(bits) != 3:
        raise Exception("Spec %s must either be a gene name, or of the form gene:pdbid:chain")
    return [bits[0].upper()] + bits[1:]

genespecs = map(parse_genespec, sys.argv[3:])
if (not have_pdbtable) and any([len(spec) == 1 for spec in genespecs]):
    raise Exception("Must specify explicit PDB IDs and chains in manual mode")
elif have_pdbtable and any([len(spec) == 3 for spec in genespecs]):
    raise Exception("Can only specify explicit PDB IDs in manual mode")

# Record a mapping from gene-name -> pdbid/chain if we're not using the database to provide that information.

if not have_pdbtable:

    genespecs_dict = dict()

    for spec in genespecs:
        if spec[0] not in genespecs_dict:
            genespecs_dict[spec[0]] = []
        genespecs_dict[spec[0]].append(spec[1:])

# Retrieve the gene records the user asked for, or all records with a PDB model but no current annotations.

if have_pdbtable:
    filters = [("%s.pdbid is not null" % pdbtable), ('%s.pdbid != "NOTFOUND"' % pdbtable)]
else:
    filters = []

if len(sys.argv) > 3:
    filters.append("(%s)" % " or ".join(["gene = '%s'" % x for x in [spec[0] for spec in genespecs]]))
else:
    # A curious arbitrary two motifs to check for. Can't remember why these two.
    filters.append("(dfg_motif_offset is null or dxxxxg_motif_offset is null)")

query = 'select gene, %s.pdbid, %s.pdbchain, offset_canonical from genes join kinases on genes.id = gene_id join inactive_pdbs on kinases.id = inactive_pdbs.id where %s' % (pdbtable, pdbtable, ' and '.join(filters))

db.query(query)

for gene, dbpdbid, dbpdbchain, offset in db.store_result().fetch_row(maxrows=0):
    
    if have_pdbtable:
        pdbs_and_chains = [(dbpdbid, dbpdbchain)]
    else:
        pdbs_and_chains = genespecs_dict[gene.upper()]

    for pdbid, pdbchain in pdbs_and_chains:

        # The model file is used to find motifs rather than relying on the sequence in the database because there are a surprising
        # number of disagreements between sequences given in Uniprot and those in PDB model files, and the model will be the authoritative
        # source of information when constructing the initial protein structure for simulation.

        # Run find-kinase-features.py to try to identify features in this kinase.

        pdbfile = os.path.join(workdir, pdbid)
        if not os.path.exists(pdbfile):
            wget_cmd = ["/usr/bin/wget", "-O", pdbfile, "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=%s" % pdbid]
            try:
                subprocess.check_call(wget_cmd, stdout=devnull, stderr=devnull)
            except Exception as e:
                raise Exception("Fetch failed: command was " + " ".join(wget_cmd))
        find_features_cmd = [os.path.join(mydir, "find-kinase-features.py"), pdbfile, pdbchain, gene, str(int(offset) + 1), workfile]
        subprocess.call(find_features_cmd)#, stdout=devnull, stderr=devnull)

        try:
            with open(workfile, "r") as f:
                results = json.load(f)
        except Exception as e:
            raise Exception("Failed to decode json resulting from %s" % " ".join(find_features_cmd))

        def offset_to_aastr(off):
            if off["match"] is not None:
                return str(off["match"] + 1)
            else:
                return off["exn"]
        def get_seq(result):
            if "seq" in result:
                return result["seq"]
            else:
                return "None"

        rec = {"Gene": gene,
               "PDB ID": pdbid,
               "Kinase offset": str(offset),
               "APE Start": offset_to_aastr(results["ape"]),
               "DFG Start": offset_to_aastr(results["dfg"]),
               "HRD Start": offset_to_aastr(results["hrd"]),
               "GxGxxG Start": offset_to_aastr(results["gxgxxg"]),
               "VAIK Start": offset_to_aastr(results["vaik"]),
               "VAIK Sequence": get_seq(results["vaik"]),
               "Activation Loop Start": str(results["activation_loop"][0]) if results["activation_loop"] is not None else "Not found",
               "Activation Loop End": str(results["activation_loop"][1]) if results["activation_loop"] is not None else "Not found",
               "aC Helix Start": str(results["ac_hel"]["start"]) if results["ac_hel"] is not None else "Not found",
               "aC Helix End": str(results["ac_hel"]["end"]) if results["ac_hel"] is not None else "Not found",
               "Bridge Glutamine": str(results["closest_glut"]),
               "Bridge Lysine": str(results["bridge_end"]),
               "Bridge Closest Atoms": str(results["closest_glut_dist"]),
               "DxxxxG Start": offset_to_aastr(results["dxxxxg"]),
               "R-Spine B4": offset_to_aastr(results["rspine_start"]),
               "R-Spine aC": offset_to_aastr(results["rspine_ac"])
        }

        writer.writerow(rec)

if workdirtemp:
    shutil.rmtree(workdir)
