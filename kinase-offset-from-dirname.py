#!/usr/bin/env python

import _mysql
import sys
import os.path

if len(sys.argv) != 2:
    print >>sys.stderr, "Usage: kinase-offset-from-dirname.py dirname"
    sys.exit(1)

base = os.path.basename(sys.argv[1])

bits = base.split("-")
if len(bits) == 2:

    print bits[1]

elif len(bits) == 1:

    db = _mysql.connect(host = "acbbdb1", user = "nsmd", passwd = "nsmdP123", db = "nsmd")
    db.query("select offset_canonical from genes join kinases on gene_id = genes.id where gene = '%s'" % base)
    rows = db.store_result().fetch_row(maxrows=0)

    if len(rows) == 0:
        raise Exception("No such gene %s" % base)
    elif len(rows) > 1:
        raise Exception("Unexpected response %s" % rows)

    print rows[0][0]

else:

    raise Exception("Can't parse %s as a gene/offset description" % base)
