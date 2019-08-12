#!/usr/bin/python

import sys
import _mysql
import difflib

db = _mysql.connect(host = "acbbdb1.picr.man.ac.uk", user = "nsmd", passwd = "nsmdP123", db = "nsmd")

if len(sys.argv) < 2:
    print >>sys.stderr, "Usage: delete-genes.py gene_name [gene_name ...]"
    sys.exit(1)

gene_spec = " or ".join(["gene = '%s'" % x for x in sys.argv[1:]])

check_q = db.query("select gene from genes inner join kinases on genes.id = gene_id inner join sequences on kinaseid = kinases.id where model_status != 0 and (%s)" % gene_spec)
rows = db.store_result().fetch_row(maxrows=0)
if len(rows) != 0:

    print >>sys.stderr, "Won't delete genes with sequences already submitted:"
    for row in rows:
        print >>sys.stderr, row[0]
    sys.exit(1)

delete_q = db.query("delete from genes, kinases, sequences using genes inner join kinases on genes.id = gene_id inner join sequences on kinases.id = kinaseid where %s" % gene_spec)
db.store_result()
