#!/usr/bin/python

import sys
import _mysql
import difflib
import importsequences
import uniprotquery
import xml.etree.ElementTree as ET

if len(sys.argv) < 2:
	print >>sys.stderr, "Usage: import-uniprot-sequences.py genesfile.csv"
	sys.exit(1)

db = _mysql.connect(host = "acbbdb1.picr.man.ac.uk", user = "nsmd", passwd = "nsmdP123", db = "nsmd")

with open(sys.argv[1], "r") as f:

        for l in f:

                bits = [x.strip() for x in l.split(",")]
                gene = bits[0]
                uniprotid = bits[1]

                fastadoc = uniprotquery.get("www.uniprot.org", "/uniprot/%s.fasta" % uniprotid)
                fastalines = [x.strip() for x in fastadoc.split("\n")]
                sequence = "".join(filter(lambda x: len(x) > 0 and x[0] != '>', fastalines))

                xmldoc = uniprotquery.get("www.uniprot.org", "/uniprot/%s.xml" % uniprotid)
                xml = ET.fromstring(xmldoc)

                kinases = []

                entry = xml.find("{http://uniprot.org/uniprot}entry")
                for feature in entry.findall("{http://uniprot.org/uniprot}feature"):

                        if feature.attrib["type"] == "domain" and feature.attrib["description"].startswith("Protein kinase"):
                                
                                loctag = feature.find("{http://uniprot.org/uniprot}location")
                                begintag = loctag.find("{http://uniprot.org/uniprot}begin")
                                endtag = loctag.find("{http://uniprot.org/uniprot}end")

                                # 1-based, inclusive indexing
                                begin = int(begintag.attrib["position"]) - 1
                                end = int(endtag.attrib["position"])
                                kinases.append(sequence[begin:end])

                if len(kinases) == 0:

                        print >>sys.stderr, "No kinases found for", gene
                        continue

                elif len(kinases) > 1:

                        print >>sys.stderr, "Multiple kinases domains registered for", gene
                        
                importsequences.import_gene(db, gene, sequence, kinases)
