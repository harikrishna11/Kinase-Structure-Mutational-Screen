#!/usr/bin/env python

import csv
import sys
import _mysql

dry_run = False

if len(sys.argv) < 3:
    print >>sys.stderr, "Usage: import-kinase-feature-spreadsheet.py kinases.csv auto|manual [force]"
    print >>sys.stderr, "auto or manual annotates whether these results were automatically generated or manually figured out, as a debugging aid. With 'force' existing records will be overwritten."
    sys.exit(1)

header_map = {
    "Gene": "gene",
    "Kinase offset": "offset",
    "APE Start": "ape_motif_offset",
    "DFG Start": "dfg_motif_offset",
    "HRD Start": "hrd_motif_offset",
    "GxGxxG Start": "gxgxxg_motif_offset",
    "VAIK Start": "vaik_motif_offset",
    "Activation Loop Start": "activation_loop_start",
    "Activation Loop End": "activation_loop_end",
    "aC Helix Start": "achelix_start",
    "aC Helix End": "achelix_end",
    "Bridge Glutamine": "bridge_glut_offset",
    "Bridge Lysine": "bridge_lys_offset",
    "Bridge Closest Atoms": "bridge_distance",
    "DxxxxG Start": "dxxxxg_motif_offset",
    "R-Spine B4": "rspine_start_offset",
    "R-Spine aC": "rspine_ac_offset"
}

query_attributes = list(header_map.values())
query_attributes.remove("gene")
query_attributes.remove("offset")

db = _mysql.connect(host = "acbbdb1.picr.man.ac.uk", user = "nsmd", passwd = "nsmdP123", db = "nsmd")

if sys.argv[2] == "auto":
    manual = False
elif sys.argv[2] == "manual":
    manual = True
else:
    raise Exception("Second argument: Must specify either auto or manual")

force = False
if len(sys.argv) >= 4 and sys.argv[3] == "force":
    force = True

with open(sys.argv[1], "r") as f:
    sheet_lines = [l.strip() for l in f.readlines()]

if len(sheet_lines) < 2:
    raise Exception("There are less than two lines in the given CSV. It should have a header (Gene,Kinase offset,...) and at least one data line")

def translate_dict(d):
    ret = dict()
    if "Gene" not in d:
        raise Exception("The given CSV file doesn't have a 'Gene' column. Is the header line missing?")
    for sheet_key, db_key in header_map.iteritems():
        if sheet_key in d:
            try:
                if db_key != "gene":
                    test = float(d[sheet_key])
                ret[db_key] = d[sheet_key]
            except ValueError:
                pass
    return ret

check_seqs = [("DFG", "dfg_motif"), ("APE", "ape_motif"), ("HRD", "hrd_motif"), ("GxGxxG", "gxgxxg_motif"), ("VAIK", "vaik_motif"), ("Bridge-glut", "bridge_glut", 1), ("Bridge-lys", "bridge_lys", 1), ("DxxxxG", "dxxxxg_motif")]

def get_check_query(s):
    motif_len = s[2] if len(s) == 3 else len(s[0])
    return "mid(canonical_sequence, %s_offset, %d)" % (s[1], motif_len)

check_queries = ", ".join(map(get_check_query, check_seqs))

def db_get(gene, offset):

    db.query("select " + ", ".join(query_attributes) + " from genes join kinases on gene_id = genes.id where gene = '%s' and offset_canonical = %d" % (gene, offset))
    rows = db.store_result().fetch_row(maxrows=0)
    if len(rows) == 0 or len(rows) > 1:
        raise Exception("Expected one record for %s / %d" % (gene, offset))

    row = rows[0]
    return dict(zip(query_attributes, row))

def db_get_sequences(gene, offset):

    db.query("select %s from genes join kinases on gene_id = genes.id where gene = '%s' and offset_canonical = %d" % (check_queries, gene, offset))
    rows = db.store_result().fetch_row(maxrows=0)
    if len(rows) == 0 or len(rows) > 1:
        raise Exception("Expected one record for %s / %d" % (gene, offset))

    row = rows[0]
    return zip([s[0] for s in check_seqs], row)

def db_update(gene, offset, attributes):

    if len(attributes) == 0:
        return

    newattribs = []
    for k, v in attributes:
        if k == "bridge_distance":
            continue
        bits = k.split("_")
        bits = bits[:-1] + ["manual"]
        if bits not in [x[0] for x in newattribs]:
            newattribs.append(("_".join(bits), "1" if manual else "0"))

    attributes.extend(newattribs)

    set_stmt = ", ".join(["%s = '%s'" % (k, v) for (k, v) in attributes])

    qry = "update genes join kinases on gene_id = genes.id set " + set_stmt + " where gene = '%s' and offset_canonical = %d" % (gene, offset)
    if dry_run:
        print >>sys.stderr, qry
    else:
        db.query(qry)
        db.store_result()

for l in csv.DictReader(sheet_lines, delimiter = ","):

    file_attribs = translate_dict(l)
    gene, offset = file_attribs["gene"], int(file_attribs["offset"])
    db_attribs = db_get(gene, offset)

    filtered_attribs = []
    
    for (k, v) in file_attribs.iteritems():
        
        if k == "gene" or k == "offset":
            continue

        dbval = db_attribs[k]
        if dbval is None or force:
            filtered_attribs.append((k, v))

    if len(filtered_attribs) == 0:
        continue

    db_update(file_attribs["gene"], int(file_attribs["offset"]), filtered_attribs)

    if manual:
        print >>sys.stderr, "Cross-check for %s/%d: annotated sequences are:" % (gene, offset)
        seqs = db_get_sequences(gene, offset)
        for (seqname, seq) in seqs:
            print >>sys.stderr, " %s: %s" % (seqname, seq)

 
