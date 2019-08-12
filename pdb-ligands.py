
import sys
import pdbquery

if len(sys.argv) < 2:
	print >>sys.stderr, "Usage: pdb-ligands.py gene-name"
	sys.exit(1)

print "\t".join(["gene", "accession_code", "resolution", "ligands"])

hits = pdbquery.getpdbinfo(sys.argv[1])

for code, desc in hits.iteritems():
    for lig in desc["ligands"]:
        print "\t".join([sys.argv[1], code, desc["resolution"], lig])
