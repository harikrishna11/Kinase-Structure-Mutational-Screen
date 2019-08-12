#!/usr/bin/python

import sys
import _mysql
import uniprotquery
import csv

db = _mysql.connect(host = "acbbdb1.picr.man.ac.uk", user = "nsmd", passwd = "nsmdP123", db = "nsmd")
db.query("select gene, kinases.id, given_sequence, matched_sequence, canonical_sequence, offset_given, offset_matched, offset_canonical, length, matched_isoform from genes inner join kinases on genes.id = kinases.gene_id where offset_canonical is not null")

for gene, kinaseid, givenseq, matchedseq, cseq, offgiven, offmatched, offc, klen, isoform in db.store_result().fetch_row(maxrows=0):

        klen = int(klen)
        offgiven = int(offgiven)
        offmatched = int(offmatched)
        offc = int(offc)

        domgiven = givenseq[offgiven : offgiven + klen]
        dommatched = matchedseq[offmatched : offmatched + klen]
        domc = cseq[offc : offc + klen]

        if domgiven != dommatched:
                print gene, kinaseid, "differs from kinase domain in closest matching Uniprot match %s" % isoform
        if dommatched != domc:
                print gene, kinaseid, "matched kinase domain differs from canonical version"


