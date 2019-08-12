#!/usr/bin/python

import sys
import httplib
import json
import copy

if len(sys.argv) < 2:
	print >>sys.stderr, "Usage: ensembl-somatic.py gene-name"
	sys.exit(1)

def get(server, relurl):

	conn = httplib.HTTPConnection(server)
	conn.request("GET", relurl)
	r1 = conn.getresponse()

	if r1.status != 200:
		print >>sys.stderr, "Failed to retrieve", "http://%s%s" % (server, relurl), "(status %d)" % r1.status
		print >>sys.stderr, r1.read()
		sys.exit(1)

	js = r1.read()
	return json.loads(js)

gene = get("beta.rest.ensembl.org", "/lookup/symbol/human/%s?content-type=application/json" % sys.argv[1])

region = "%s:%s-%s" % (gene["seq_region_name"], gene["start"], gene["end"])

print >>sys.stderr, "Looking up", gene["description"], "at", region

mutations = get("beta.rest.ensembl.org", "/feature/region/human/%s?feature=somatic_variation&content-type=application/json" % region)

mutation_consequences = []

for mut in mutations:

	effects = get("beta.rest.ensembl.org", "/vep/human/id/%s/consequences?content-type=application/json" % mut["ID"])

	for trans in effects["data"][0]["transcripts"]:

		if trans["biotype"] != "protein_coding":
			continue

		if trans["translation_start"] is not None:
			
			ts = trans["translation_start"]
			te = trans["translation_end"]

			if ts != te:
				this_aa = "%d-%d" % (ts, te)
			else:
				this_aa = str(ts)

			for a in trans["alleles"]:

				pep_var = a["pep_allele_string"]
				if pep_var is not None:

					mut_con = copy.deepcopy(mut)
					mut_con["pep_alleles"] = pep_var
					mut_con["aa_index"] = this_aa
					mut_con["consequence_type"] = str(a["consequence_terms"])
					mut_con["alt_alleles"] = a["allele_string"]
					mut_con["transcript_id"] = trans["transcript_id"]
					mut_con["gene_id"] = trans["gene_id"]
		
					mutation_consequences.append(mut_con)

seen_keys = set()

for mut in mutation_consequences:
	for k, v in mut.iteritems():
		if k not in seen_keys:
			seen_keys.add(k)

seen_keys = list(seen_keys)

print "\t".join(seen_keys)

for mut in mutation_consequences:
	print "\t".join([str(mut[k]) for k in seen_keys])
