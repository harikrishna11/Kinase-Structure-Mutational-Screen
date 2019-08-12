#!/apps/scripts/nsmd/deps/modeller-9.14-install/bin/modpy.sh python

import modeller
import modelutils
import sys
import re
import json
import _mysql

## This script tries to automatically find interesting kinase features in a given PDB file / chain.
## A gene and kinase domain start coordinate are also needed for cross-checking with the database,
## and to use existing manual (user-supplied) motifs when available.
## It can either output user-friendly information or a JSON description for consumption by report-all-kinase-features.py.

if len(sys.argv) < 5:
    print >>sys.stderr, "Usage: find-kinase-features query.pdb chain gene kinase_domain_start [json]"
    sys.exit(1)

db = _mysql.connect(host = "acbbdb1.picr.man.ac.uk", user = "nsmd", passwd = "nsmdP123", db = "nsmd")

pdbfile = sys.argv[1]
chain = sys.argv[2]
gene = sys.argv[3]
kinase_domain_start = int(sys.argv[4])
if len(sys.argv) >= 6:
    json_out = sys.argv[5]
else:
    json_out = None

with open("/dev/null", "w") as fnull:

    # Supress verbose version notice
    oldout = sys.stdout
    sys.stdout = fnull
    env = modeller.environ()
    mdl = modeller.model(env, file = pdbfile, model_segment = ('FIRST:' + chain, 'LAST:' + chain))
    sys.stdout = oldout

# Read the feature annotations given in the PDB file.

seq = modelutils.model_to_sequence(mdl)

hels = modelutils.read_helices(pdbfile, chain)
sheets = modelutils.read_sheets(pdbfile, chain)

# Fetch existing motif offsets if known, to find derived features.
def intornone(x):
    if x is None:
        return x
    return int(x)

db.query("select ape_motif_offset, dfg_motif_offset, vaik_motif_offset from genes join kinases on gene_id = genes.id where gene = '%s' and offset_canonical = %d" % (gene, kinase_domain_start - 1))
oldape, olddfg, oldvaik = map(intornone, db.store_result().fetch_row(maxrows=0)[0])

# Check whether the PDB file's given relationship to a Uniprot ID is consistent, by comparing DBREF
# lines in the PDB against a simple sequence alignment.

pdb_to_uniprot = modelutils.read_pdb_to_uniprot(pdbfile, chain)
aligned_uniprot_offset = modelutils.align_pdb_to_uniprot(gene, kinase_domain_start, seq, db)

# Find the most common offset in the sequence alignment.
# Positive aligned_uniprot_offset means the PDB file is ahead of the Uniprot coordinates.

count_offsets = dict()
for (pdboff, uniprotoff) in pdb_to_uniprot.iteritems():
    diff = pdboff - uniprotoff
    if diff not in count_offsets:
        count_offsets[diff] = 0
    count_offsets[diff] += 1

most_common_offset = max(count_offsets.iteritems(), key = lambda x : x[1])[0]

# It is distressingly common for PDB files' described Uniprot mappings to be off-by-one,
# or described with respect to a non-canonical isoform.

if most_common_offset != aligned_uniprot_offset:
    print >>sys.stderr, "Sequence alignment and DBREF lines disagree (%d vs. %d)" % (aligned_uniprot_offset, most_common_offset)
    print >>sys.stderr, "Assuming PDB error; adjusting offset dict. Check results carefully."

    for k in pdb_to_uniprot:
        pdb_to_uniprot[k] += (most_common_offset - aligned_uniprot_offset)

uniprot_to_pdb = dict([(v, k) for (k, v) in pdb_to_uniprot.iteritems()])

kinase_domain_start = uniprot_to_pdb[kinase_domain_start]

# OK now to find the interesting motifs:

# Find a unique regex match or verbosely report a failure to do so:
def findone(motif, motifname, offset = None, limit = None):

    if offset is not None:
        if limit is not None:
            tosearch = seq[offset:limit]
        else:
            tosearch = seq[offset:]
    else:
        tosearch = seq
        offset = 0
    matches = list(re.finditer(motif, tosearch))
    if len(matches) == 0:
        raise Exception("Failed to find " + motifname + " motif")
    elif len(matches) >= 2:
        raise Exception("Motif " + motifname + " is ambiguous: matches start at " + ", ".join([str(m.start() + offset + 1) for m in matches]))
    return matches[0]

# Report a motif if found. The "use old version from the database" logic applies when the user has manually
# given e.g. the location of the correct APE motif, which we can then use to auto-find derived structures
# (e.g. the activation loop, whose location is defined relative to other motifs)
def report_motif(motif, motifname = "", startoff = None, limitoff = None, old = None):

    if motifname == "":
        motifname = motif
    if old is not None:
        old = uniprot_to_pdb[old - 1]
    try:
        match = findone(motif, motifname, offset = startoff, limit = limitoff)
        if startoff is None:
            startoff = 0
        if old is not None and old != (match.start() + startoff):
            raise Exception("Search for %s produced %d which conflicts with existing finding %d" % (motifname, match.start() + startoff, old))
    except Exception as e:
        if old is not None:
            print >>sys.stderr, e
            print >>sys.stderr, "Using old value from database: %s / %d" % (motifname, old)
            return {"match": old, "seq": None, "exn": None}
        else:
            return {"match": None, "exn": str(e)}
    return {"match": match.start() + startoff, "seq": match.group(0), "exn": None}

# Find each simple moif we're interested in:

ape = report_motif("APE", old = oldape)
dfg = report_motif("DFG", old = olddfg)
hrd = report_motif("HRD")
gxgxxg = report_motif("G.G..G", motifname = "GxGxxG")
vaik = report_motif("VA[IV]K", motifname = "VAIK/VAVK", old = oldvaik)
if vaik["exn"] is not None:
    vaik2 = report_motif(".A.K", motifname = ".A.K")
    if vaik2["exn"] is not None:
        vaik2["exn"] = "%s; %s" % (vaik["exn"], vaik2["exn"])
    vaik = vaik2

# Careful of 1-based helix reports and 0-based string offsets here.
if ape["match"] is not None:
    helices_after_ape = [h for h in hels if (h["start"] - 1) > ape["match"] and h["class"] == 1]
    if len(helices_after_ape) != 0:
        limit = helices_after_ape[0]["end"] - 1
        dxxxxg = report_motif("D....G", motifname = "DxxxxG", startoff = ape["match"], limitoff = limit)
    else:
        dxxxxg = {"match": None, "exn": "No helices found beyond APE offset " + ape["match"]}
else:
    dxxxxg = {"match": None, "exn": "Can't find DxxxxG without APE"}

activation_loop = None

# dfg is 0-based, hel is 1-based
if dfg["match"] is not None:
    for hel in hels:
        if hel["start"] >= dfg["match"] + 4 and hel["class"] == 1:
            activation_loop = (dfg["match"] + 4, hel["start"])
            break

# Find the 'K' in VAI/VK
if vaik["match"] is not None:
    bridge_end = [res for res in mdl.residues if int(res.num) == vaik["match"] + 3 + 1][0] # + 1 for 1-based coordinates
else:
    bridge_end = None

# Ad-hoc measure of helix 'niceness', for finding the alpha-C helix. We're looking for one that
# contains an E residue that's close to the 'K' of the VAIK motif.
def score_helix(hel):

    gluts_in_hel = [(res, modelutils.closest_res_dist(res, bridge_end)) for res in mdl.residues if int(res.num) >= hel["start"] and int(res.num) <= hel["end"] and res.code == "E"]
    if len(gluts_in_hel) == 0:
	return (None, 9999)
    closest_glut, closest_glut_dist = min(gluts_in_hel, key = lambda x: x[1])
    return (closest_glut, closest_glut_dist)

if bridge_end is not None:

    long_hels_in_kinase = [hel for hel in hels if hel["start"] >= kinase_domain_start and int(hel["end"]) - int(hel["start"]) >= 5]
    scored_hels = sorted([(hel, score_helix(hel)) for hel in long_hels_in_kinase], key = lambda x : x[1][1])
    (ac_hel, (closest_glut, closest_glut_dist)) = scored_hels[0]
    search_offset = ac_hel["end"] - 2
    search_limit = search_offset + 3
    rspine_ac = report_motif("[LIVM]", motifname = "aC helix R-spine", startoff = search_offset, limitoff = search_limit)

else:
        
    ac_hel, closest_glut, closest_glut_dist = None, None, None
    rspine_ac = {"match": None, "exn": "Need bridge end to search for the aC helix"}

sheets = [s for s in sheets if s["start"] >= kinase_domain_start]

# Find part of the R-spine. This should be in the 4th beta-sheet, but this is pretty unreliable
# and user checking is highly recommended.

if len(sheets) >= 4:
    sheetstart = sheets[3]["start"]
    search_offset = sheetstart
    search_limit = search_offset + 1
    rspine_start = report_motif("[LIVFY]", motifname = "Beta-4 R-spine", startoff = search_offset, limitoff = search_limit)
else:
    rspine_start = {"match": None, "exn": "Not enough beta-sheet records in kinase domain to find R-spine start"}

# Functions to translate PDB coordinates into Uniprot ones.
def to_up(idx):
    try:
        return pdb_to_uniprot[int(idx)]
    except Exception as e:
        print >>sys.stderr, "Warning: could not convert", idx, e
        return idx

def motif_to_up(m):

    if m is None:
        return None
    m["match"] = to_up(m["match"])
    return m

# Translate everything from PDB to Uniprot coordinates before output
ape = motif_to_up(ape)
dfg = motif_to_up(dfg)
hrd = motif_to_up(hrd)
gxgxxg = motif_to_up(gxgxxg)
dxxxxg = motif_to_up(dxxxxg)
vaik = motif_to_up(vaik)
rspine_start = motif_to_up(rspine_start)
rspine_ac = motif_to_up(rspine_ac)
if activation_loop is not None:
    activation_loop = (to_up(activation_loop[0]), to_up(activation_loop[1]))
if bridge_end is not None:
    bridge_end = to_up(bridge_end.num)
if ac_hel is not None:
    ac_hel["start"] = to_up(ac_hel["start"])
    ac_hel["end"] = to_up(ac_hel["end"])
if closest_glut is not None:
    closest_glut = to_up(closest_glut.num)

if json_out is not None:

    with open(json_out, "w") as f:
        report = {"ape": ape, 
                  "dfg": dfg, 
                  "hrd": hrd, 
                  "gxgxxg": gxgxxg, 
                  "dxxxxg": dxxxxg, 
                  "vaik": vaik,
                  "rspine_start": rspine_start,
                  "rspine_ac": rspine_ac,
                  "activation_loop": activation_loop, 
                  "bridge_end": bridge_end if bridge_end is not None else None, 
                  "ac_hel": ac_hel, 
                  "closest_glut": closest_glut if closest_glut is not None else None, 
                  "closest_glut_dist": closest_glut_dist}
        json.dump(report, f)

else:

    def print_match(m, mname):
        if m["match"] is not None:
            print "Motif", mname, "found starting at position", m["match"] + 1
        else:
            print "*** Motif", mname, "NOT FOUND (%s)" % m["exn"]
            
    for (m, mname) in [(ape, "APE"), (dfg, "DFG"), (hrd, "HRD"), (gxgxxg, "GxGxxG"), (vaik, "VAIK"), (dxxxxg, "DxxxxG"), (rspine_ac, "R-Spine aC"), (rspine_start, "R-Spine start")]:
        print_match(m, mname)

    if activation_loop is not None:
        print "Activation loop spans from residue", activation_loop[0], "to", activation_loop[1]

    if ac_hel is not None:
        print "aC-helix spans from residue", ac_hel["start"], "to", ac_hel["end"]
        print "Closest E in the helix to K of VAIK/VAVK: residue", closest_glut, "with a distance of", closest_glut_dist
    else:
        print "aC helix not found, or VAIK's 'K' missing"

