#!/usr/bin/python

import sys
import mutate

try:

    if len(sys.argv) < 2:
        raise
	
    old_base, new_base, base_pos = mutate.parsespec(sys.argv[1])

except:

    print >>sys.stderr, "Usage: mutate-protein.py mutation-descriptor"
    print >>sys.stderr, "Mutation descriptor should be like AnB where n is a one-based index."
    sys.exit(1)

in_bases = ""
header = None

longest_col = 0

for line in sys.stdin:

    line = line.strip()
        
    if line.startswith(">"):
        if header is not None:
            raise Exception("Saw more than one FASTA header line")
        header = line
        continue

    if len(line) > longest_col:
        longest_col = len(line)

    in_bases = in_bases + line

out_bases = mutate.applymutation(in_bases, old_base, new_base, base_pos)

if header is not None:
	sys.stdout.write(header)
	sys.stdout.write(" (mutated %s)\n" % sys.argv[1])

while len(out_bases) > 0:

	sys.stdout.write("%s\n" % out_bases[:longest_col])
	out_bases = out_bases[longest_col:]
