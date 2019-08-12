#!/usr/bin/env python

import csv
import httplib
import _mysql
import sys
import math

db = _mysql.connect(host = "acbbdb1", user = "nsmd", passwd = "nsmdP123", db = "nsmd")

db.query("select sequences.id, gene, sourceid from sequences join kinases on kinaseid = kinases.id join genes on gene_id = genes.id where sourceid != 'Wild type' and pdbid != 'NOTFOUND' and pdbid is not null")
rows = db.store_result().fetch_row(maxrows=0)

fetch_genes = []

for (idx, gene, sourceids) in rows:

    sources = filter(lambda source: source.find("LUNG") != -1, sourceids.split(","))
    if len(sources) == 0:
        continue

    if gene not in fetch_genes:
        fetch_genes.append(gene.upper())

gene_levels = dict()

with open(sys.argv[1], "r") as f:

    lines = f.readlines()

    for l in csv.DictReader(lines, delimiter = "\t"):

        gene = l["Description"].strip()
        if len(gene) == 0 or gene.isdigit():
            continue

        if gene.upper() not in fetch_genes:
            continue
        
        for (sample, level) in l.iteritems():

            if sample is None or len(sample.strip()) == 0 or sample == "Description" or sample == "Accession":
                continue

            flevel = float(level)
            if math.isnan(flevel):
                print >>sys.stderr, "Warning: NaN for %s, %s" % (gene, sample)
            gene_levels[(gene.upper(), sample)] = flevel

for (idx, gene, sourceids) in rows:

    sources = filter(lambda source: source.find("LUNG") != -1, sourceids.split(","))
    if len(sources) == 0:
        continue

    def try_get_level(gene, sample):
        try:
            ret = gene_levels[(gene, sample)]
            if math.isnan(ret):
                ret = 0.0
            return ret
        except KeyError:
            print >>sys.stderr, "No record for %s" % ((gene, sample),)
            return 0.0

    best_expr = max([try_get_level(gene.upper(), sample) for sample in sources])

    query = "update sequences set expression_level = %g where id = %d" % (best_expr, int(idx))
    
    db.query(query)
    db.store_result()
