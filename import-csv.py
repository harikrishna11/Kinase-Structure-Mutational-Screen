#!/usr/bin/python

import sys
import _mysql
import difflib
import importsequences

if len(sys.argv) < 2:
	print >>sys.stderr, "Usage: map-genes-to-uniprot.py genes.csv"
	sys.exit(1)

db = _mysql.connect(host = "acbbdb1.picr.man.ac.uk", user = "nsmd", passwd = "nsmdP123", db = "nsmd")

with open(sys.argv[1], "r") as gf:

        for l in gf:

                bits = l.split(",")
                bits = map(lambda x: x.strip(), bits)

                if len(bits) < 3:
                        print >>sys.stderr, "Parse error at", l
                        sys.exit(1)

                importsequences.import_gene(db, bits[0], bits[1], bits[2:4])

