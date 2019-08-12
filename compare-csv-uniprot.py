#!/usr/bin/python

import sys
import _mysql
import csv

db = _mysql.connect(host = "acbbdb1.picr.man.ac.uk", user = "nsmd", passwd = "nsmdP123", db = "nsmd")

def getuniprotid(gene):

        db.query("select uniprotid from genes where LOWER(gene) = '%s'" % gene.lower())
        rows = db.store_result().fetch_row(maxrows=0)
        if len(rows) == 1:
                return rows[0][0]
        elif len(rows) > 1:
                raise Exception("Multiple uniprot IDs for %s" % gene)
        else:
                raise Exception("No match")

with open(sys.argv[1], "r") as f:

        for l in f:

                bits = [x.strip() for x in l.split(",")]
                try:
                        existing_id = getuniprotid(bits[0])
                except:
                        print bits[0], "not in DB"
                        continue

                if existing_id is None:
                        print bits[0], "has no existing match, file suggests", bits[1]
                        continue

                if existing_id != bits[1]:
                        print bits[0], "mapped to", existing_id, "but file suggests", bits[1]
                        continue

                print bits[0], "agrees with DB"

with open(sys.argv[1], "r") as f:
        for l in f:
                bits = [x.strip() for x in l.split(",")]

                db.query("select gene from genes where uniprotid = '%s'" % bits[1])
                rows = db.store_result().fetch_row(maxrows=0)
                for row in rows:
                        if row[0].upper() != bits[0].upper():
                                print "Uniprot ID", bits[1], "recommended for", bits[0], "but already associated with", row[0]



