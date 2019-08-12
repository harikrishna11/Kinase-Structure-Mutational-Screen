#!/usr/bin/python

import sys
import _mysql
import uniprotquery

def checksequence(sequence, uniprotid, geneid, gene, db, fix=False):

    fastadoc = uniprotquery.get("www.uniprot.org", "/uniprot/%s.fasta" % uniprotid)
    fastalines = [x.strip() for x in fastadoc.split("\n")]
    realsequence = "".join(filter(lambda x: len(x) > 0 and x[0] != '>', fastalines))

    if sequence != realsequence:

        print "Gene", gene, "sequence does not match"

    else:
        
        print "Gene", gene, "matches"

if __name__ == "__main__":
    
    db = _mysql.connect(host = "acbbdb1.picr.man.ac.uk", user = "nsmd", passwd = "nsmdP123", db = "nsmd")
    db.query("select canonical_sequence, uniprotid, id, gene from genes where canonical_sequence is not null")
    rows = db.store_result().fetch_row(maxrows=0)
    
    fix = len(sys.argv) >= 2 and sys.argv[1] == "--fix"

    for row in rows:

        checksequence(*row, db=db, fix=fix)
