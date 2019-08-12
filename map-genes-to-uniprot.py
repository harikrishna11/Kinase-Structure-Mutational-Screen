#!/usr/bin/python

import sys
import _mysql
import uniprotquery
import csv

db = _mysql.connect(host = "acbbdb1.picr.man.ac.uk", user = "nsmd", passwd = "nsmdP123", db = "nsmd")

def getuniprotid(gene):

        db.query("select uniprotid from genes where gene = '%s'" % gene)
        rows = db.store_result().fetch_row(maxrows=0)
        if len(rows) == 1:
                return rows[0][0]
        elif len(rows) > 1:
                raise Exception("Multiple uniprot IDs for %s" % gene)
        else:
                return None

def storeuniprotid(recordid, result):

        db.query("update genes set uniprotid = '%s', matched_isoform = '%s', given_matched_score = %g, given_matched_diffops = '%s', matched_sequence = '%s' where id = %d" % (result["uniprotid"], result["isoform"], result["score"], result["diffops"], result["matched_sequence"], recordid))
        db.store_result()

db.query("select id, gene, given_sequence from genes")

for recordid, gene, given_seq in db.store_result().fetch_row(maxrows=0):

        uniprotid = getuniprotid(gene)
        if uniprotid is not None:
                continue
                
        uniprotid = uniprotquery.getbestmatch(gene, given_seq)
        if uniprotid is not None:
                storeuniprotid(int(recordid), uniprotid)

        if uniprotid is not None:
                print >>sys.stderr, gene, "->", uniprotid["isoform"]
        else:
                print >>sys.stderr, "No Uniprot results for", gene

