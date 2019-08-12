import math
import itertools
import sys
import os.path
import subprocess
import difflib

## Common model-manipulation routines, used by find-kinase-features.py, check-active-dfg.py and others.

# Read the HELIX records from a PDB file that relate to a particular chain.
# Note the PDB format uses fixed column offsets (!) *not* space-delimited fields.
def read_helices(pdbfile, chain):

    result = []

    with open(pdbfile, "r") as f:
        for l in f:
            if l.startswith("HELIX"):

                thischain = l[19]
                if thischain != chain:
                    continue

                newrec = dict()

                newrec["name"] = l[11:14]
                newrec["start"] = int(l[21:25])
                newrec["end"] = int(l[33:37])
                newrec["class"] = int(l[38:40])

                result.append(newrec)

    return result

# Read beta sheets from a PDB file. Similar caveats to the above.
def read_sheets(pdbfile, chain):

    result = []

    with open(pdbfile, "r") as f:
        for l in f:
            if l.startswith("SHEET"):
                
                thischain = l[21]
                if thischain != chain:
                    continue

                newrec = dict()

                newrec["start"] = int(l[22:26])
                newrec["end"] = int(l[33:37])
                newrec["id"] = "%s-%s" % (l[11:14], l[7:10])

                result.append(newrec)

    return result
        
# Convert a python diffops match tuple to a match length. Only used below.
def get_match_len(optuple):
    
    tag, i1, i2, j1, j2 = optuple
    return i2 - i1

# Get the uniprot sequence for a given gene/kinase offset from the database and check how it aligns to 'pdbseq',
# the base sequence extracted from a PDB model of the same kinase.
def align_pdb_to_uniprot(gene, offset, pdbseq, dbhandle):

    dbhandle.query("select canonical_sequence from genes join kinases on gene_id = genes.id where gene = '%s' and offset_canonical = %d" % (gene, offset - 1))
    uniprotseq = dbhandle.store_result().fetch_row(maxrows=0)
    if(len(uniprotseq) != 1):
        raise Exception("Expected one result for %s/%d; got %s" % (gene, offset, uniprotseq))
    uniprotseq = uniprotseq[0][0]

    alignment = difflib.SequenceMatcher(None, uniprotseq, pdbseq, autojunk = False).get_opcodes()

    match_or_mismatch_ops = filter(lambda optuple: optuple[0] == "equal", alignment)
    longest_match = max(match_or_mismatch_ops, key = get_match_len)
    return longest_match[3] - longest_match[1]        

# Find out how a PDB file asserts it aligns to a Uniprot protein.
# This is done with DBREF records in the PDB file which give the 'other database' as UNP.
def read_pdb_to_uniprot(pdbfile, chain, lowest_pdb_offset = -5000, highest_pdb_offset = 5000):

    result = dict()
    with open(pdbfile, "r") as f:
        
        unique_offset = None
        unique_offset_valid = True

        for l in f:

            if not l.startswith("DBREF"):
                continue

            try:
                tag, pdbcode, pdbchain, pdbstart, pdbend, otherdb, otherdbid, otherdbdesc, otherstart, otherend = l.split()
                pdbstart = int(pdbstart)
                pdbend = int(pdbend)
                otherstart = int(otherstart)
                otherend = int(otherend)
            except ValueError as e:
                print >>sys.stderr, "Warning: ignore malformed DBREF line", l
                print >>sys.stderr, e
                continue

            if pdbchain != chain or otherdb != "UNP":
                continue
                
            offset = otherstart - pdbstart

            if unique_offset is None:
                unique_offset = offset
            elif unique_offset != offset:
                unique_offset_valid = False
            for i in range(pdbstart, pdbend + 1):
                result[i] = i + offset

        # Assume a straightforward mapping for other offsets if possible:
        if unique_offset_valid and unique_offset is not None:
            for i in range(lowest_pdb_offset, highest_pdb_offset):
                result[i] = i + unique_offset

    return result

# Helper functions for measuring inter-atom and inter-residue distances. The residues and atoms
# are MODELLER data structures.
def euclid_dist(a1, a2):
    return math.sqrt(math.pow(a1.x - a2.x, 2) + math.pow(a1.y - a2.y, 2) + math.pow(a1.z - a2.z, 2))

def closest_res_dist(res1, res2):
    return min([euclid_dist(a1, a2) for (a1, a2) in itertools.product(res1.atoms, res2.atoms)])

def atom_dist(mdl, res1, at1, res2, at2):
    dist = euclid_dist(mdl.residues[res1].atoms[at1], mdl.residues[res2].atoms[at2])
    print >>sys.stderr, "Distance %s:%s <-> %s:%s = %g" % (mdl.residues[res1].code, at1, mdl.residues[res2].code, at2, dist)
    return dist

# Extract the base sequence from a PDB file. Rather than parsing it directly, use MODELLER's mdl data structure.
def model_to_sequence(mdl):
    
    codes = []

    for res in mdl.residues:

        # Negative-numbered bases are hopefully always tags and other extraneous stuff!
        if int(res.num) <= 0:
            continue

        while len(codes) < int(res.num) - 1: # 1-based coordinates
            codes.append("?")
        codes.append(res.code)
 
    return "".join(codes)

# Download a PDB file if necessary, and/or supply a corrected one from $THISDIR/corrected_pdbs.
def get_pdb_file(pdbid, workdir):

    pdbfile = os.path.join(workdir, pdbid)

    mydir = os.path.dirname(os.path.realpath(__file__))

    corrected_pdbs = os.path.join(mydir, "corrected_pdbs")
    corrected_pdb_file = os.path.join(corrected_pdbs, pdbid)

    if os.path.exists(corrected_pdb_file):
        print >>sys.stderr, "Notice: using corrected PDB for", pdbid
        return corrected_pdb_file
    elif os.path.exists(pdbfile):
        return pdbfile
    
    wget_cmd = ["/usr/bin/wget", "-O", pdbfile, "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=%s" % pdbid]
    try:
        print "Fetch", pdbfile
        with open("/dev/null", "w") as devnull:
            subprocess.check_call(wget_cmd, stdout=devnull, stderr=devnull)
    except Exception as e:
        raise Exception("Fetch failed: command was " + " ".join(wget_cmd) + "; exception was: " + str(e))

    return pdbfile
