#!/usr/bin/python

import sys
import httplib
import _mysql
import pdbquery
import csv
import mutate
import uniprotquery
import xml.etree.ElementTree as ET

def get(server, relurl):

	conn = httplib.HTTPConnection(server)
	conn.request("GET", relurl)
	r1 = conn.getresponse()

	if r1.status != 200:
		print >>sys.stderr, "Failed to retrieve", "http://%s%s" % (server, relurl), "(status %d)" % r1.status
		print >>sys.stderr, r1.read()
		sys.exit(1)

	return r1.read()

db = _mysql.connect(host = "acbbdb1", user = "nsmd", passwd = "nsmdP123", db = "nsmd")

def maybe_create_sequence(kid, mutation, seq, source, sourceid):

        db.query("select id, sourceid from sequences where kinaseid = '%s' and sequence = '%s'" % (kid, seq))
        rows = db.store_result().fetch_row(maxrows=0)
        if len(rows) == 1:
                sources = rows[0][1].split(",")
                if sourceid not in sources:
                        db.query("update sequences set sourceid = '%s,%s' where id = %d" % (rows[0][1], sourceid, int(rows[0][0])))
                        db.store_result()
        elif len(rows) > 1:
                raise Exception("Duplicate sequences %s" % seq)
        else:
                db.query("insert into sequences (kinaseid, mutation, sequence, source, sourceid) values ('%s', '%s', '%s', '%s', '%s')" % (kid, mutation, seq, source, db.escape_string(sourceid)))
                db.store_result()

def get_gene(name, count = True):

        db.query("select id from genes where name = '%s'" % name)
        rows = db.store_result().fetch_row(maxrows=0)
        if len(rows) == 1:
                if count:
                        genes_existing += 1
                return int(rows[0][0])
        elif len(rows) > 1:
                raise Exception("Duplicate rows for gene %s!" % name)
        else:
                return None

db.query("select kinases.id, gene, canonical_sequence, offset_canonical, length, uniprotid from genes inner join kinases on genes.id = kinases.gene_id where offset_canonical is not null and pdbid is not null and pdbid != 'NOTFOUND'")
rows = db.store_result().fetch_row(maxrows=0)

rows_by_gene = {row[1].upper(): row for row in rows}

genelist = ",".join([row[1] for row in rows])

cbio_tsv = get("www.cbioportal.org", "/webservice.do?cmd=getMutationData&cancer_study_id=cellline_ccle_broad&genetic_profile_id=cellline_ccle_broad_mutations&gene_list=%s" % genelist)
cbio_lines = cbio_tsv.split("\n")

while cbio_lines[0].startswith("#"):
        print >>sys.stderr, "cBio search warning:", cbio_lines[0]
        cbio_lines = cbio_lines[1:]

mut_errors = dict() 

for row in rows:

        domain = row[2][int(row[3]) : int(row[3]) + int(row[4])]
        maybe_create_sequence(kid = row[0], mutation = "None", seq = domain, source = "Wild type", sourceid = "%s_%s_wild" % (row[1], row[0]))

tried_mutations = dict()

for l in csv.DictReader(cbio_lines, delimiter = "\t"):

        # Skip non-missense mutations for now: don't know how to apply cBio's syntax.
        if l["mutation_type"] != "Missense_Mutation":
                continue

        gene = l["gene_symbol"]
        try:
                kinaseid, gene, baseseq, domain_offset, kinase_length, uniprotid = rows_by_gene[gene.upper()]
                domain_offset = int(domain_offset)
                kinase_length = int(kinase_length)
        except KeyError:
                print >>sys.stderr, gene, "not found?"

        domain = baseseq[domain_offset : domain_offset + kinase_length]

        if gene not in mut_errors:
                mut_errors[gene] = {"muts": 0, "errors": []}

        try:

                oldbase, newbase, pos = mutate.parsespec(l["amino_acid_change"])

                # pos is a 1-based coordinate, hence the extra +/- 1s here.
                if (pos + len(newbase) - 2) < domain_offset or (pos - 1) >= domain_offset + len(domain):
                        continue

                if gene not in tried_mutations:
                        tried_mutations[gene] = []
                tried_mutations[gene].append(l)
        
                mut_errors[gene]["muts"] += 1        
                mutseq = mutate.applymutation(baseseq, oldbase, newbase, pos)

        except Exception as e:

                error_report = "Ignored %s mutation %s: %s" % (gene, l["amino_acid_change"], e)
                mut_errors[gene]["errors"].append(error_report)
                continue

        mutdesc = l["amino_acid_change"]

        source = "cBio %s" % l["genetic_profile_id"]
        sourceid = l["case_id"]

        maybe_create_sequence(kinaseid, mutdesc, mutseq[domain_offset : domain_offset + kinase_length], source, sourceid)

for gene, errors in mut_errors.iteritems():

        if len(errors["errors"]) == 0:
                continue

        print >>sys.stderr, "Gene %s: %d/%d mutations failed. Detail:" % (gene, len(errors["errors"]), errors["muts"])
        for e in errors["errors"]:
                print >>sys.stderr, e

        kinaseid, gene, baseseq, domain_offset, kinase_length, uniprotid = rows_by_gene[gene.upper()]
        
        uniprotdoc = uniprotquery.get("www.uniprot.org", "/uniprot/%s.xml" % uniprotid)
        upxml = ET.fromstring(uniprotdoc)
        result = upxml.find("{http://uniprot.org/uniprot}entry")
        isoforms = uniprotquery.getisoforms(result, baseseq, "Base")

        if len(isoforms) != 1:

                for isoname, isoform in isoforms.iteritems():

                        if isoform == baseseq:
                                continue

                        worked = 0

                        for l in tried_mutations[gene]:

                                oldbase, newbase, pos = mutate.parsespec(l["amino_acid_change"])
                                try:
                                        mutate.applymutation(isoform, oldbase, newbase, pos)
                                        worked += 1
                                except:
                                        pass

                        print "Note: %d/%d worked for isoform %s" % (worked, len(tried_mutations[gene]), isoname)
                                                                     
