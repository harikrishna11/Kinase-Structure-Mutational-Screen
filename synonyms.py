#!/usr/bin/python

import sys
import urllib2
import urllib

def get(server, relurl, querydoc = None):

	req = urllib2.Request("http://%s%s" % (server, relurl), data=querydoc)

	try:

		f = urllib2.urlopen(req)
		return f.read()

	except Exception as e:
		print >>sys.stderr, "Failed to retrieve", "http://%s%s" % (server, relurl), e
		print >>sys.stderr, "Query doc:"
		print >>sys.stderr, querydoc
		print >>sys.stderr, "\nResponse:"
		print >>sys.stderr, e.read()

f = open(sys.argv[1], "r")

for l in f:

        bits = l.split(",")
        gene = bits[0].strip()

        querydoc = get("www.bioinformatics.org", "/textknowledge/synonym.php", urllib.urlencode({"term": gene, "tax": "9606"}))

        bits = querydoc.split()
        headerskipped = False

        synonyms = []

        for i in range(len(bits)):

                if bits[i] == "**":
                        if not headerskipped:
                                headerskipped = True
                                continue
                        synonyms.append(bits[i-1])

        if len(synonyms) == 0:
                print "%s,???" % gene
        elif gene in synonyms:
                print "%s,%s" % (gene, gene)
        else:
                print "%s,%s" % (gene, synonyms[0])

