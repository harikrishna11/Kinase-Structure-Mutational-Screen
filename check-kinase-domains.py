#!/usr/bin/python

import sys
import _mysql
import uniprotquery
import xml.etree.ElementTree as ET

def checkdomain(offset, length, uniprotid, gene, kinaseid, db, fix=False):

        offset = int(offset)
        length = int(length)

        xmldoc = uniprotquery.get("www.uniprot.org", "/uniprot/%s.xml" % uniprotid)
        xml = ET.fromstring(xmldoc)

        best = None
        besttup = None

        entry = xml.find("{http://uniprot.org/uniprot}entry")
        for feature in entry.findall("{http://uniprot.org/uniprot}feature"):

                if feature.attrib["type"] == "domain" and (feature.attrib["description"].startswith("Protein kinase") or feature.attrib["description"].startswith("Histidine kinase")):
                                
                        loctag = feature.find("{http://uniprot.org/uniprot}location")
                        begintag = loctag.find("{http://uniprot.org/uniprot}begin")
                        endtag = loctag.find("{http://uniprot.org/uniprot}end")

                        # 1-based, inclusive indexing
                        dboffset = int(begintag.attrib["position"]) - 1
                        dblength = int(endtag.attrib["position"]) - dboffset
                        
                        diff = abs(dboffset - offset) + abs(dblength - length)
                        if best is None or best > diff:
                                best = diff
                                besttup = (dboffset, dblength)

        if best is None:
                print >>sys.stderr, "Uniprot", uniprotid, "for gene", gene, "does not mention a kinase domain"
        elif best > 0:
                db.query("update kinases set offset_canonical = %d, length = %d where id = %s" % (besttup[0], besttup[1], kinaseid))
                db.store_result()
                db.query("delete from sequences where kinaseid = %s" % kinaseid)
                db.store_result()
                print >>sys.stderr, "Uniprot", uniprotid, "for gene", gene, "error magnitude", best, "Uniprot says (%d-%d), Excel sheet said (%d-%d) %s" % (offset+1, offset+length, besttup[0] + 1, besttup[0] + besttup[1], " (fixed)" if fix else "")
        else:
                print >>sys.stderr, "Uniprot", uniprotid, "for gene", gene, "matches"

if __name__ == "__main__":

	db = _mysql.connect(host = "acbbdb1.picr.man.ac.uk", user = "nsmd", passwd = "nsmdP123", db = "nsmd")
	db.query("select offset_canonical, length, uniprotid, gene, kinases.id from genes inner join kinases on gene_id = genes.id where offset_canonical is not null and pdbid is not null and pdbid != 'NOTFOUND'")
	rows = db.store_result().fetch_row(maxrows=0)

	fix = len(sys.argv) >= 2 and sys.argv[1] == "--fix"

	for row in rows:

		checkdomain(*row, db=db, fix=fix)

