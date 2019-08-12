#!/usr/bin/env python

# For each entry in the kinase database with a model specified, check if the model disagrees with the sequence at any point.

import csv
import httplib
import _mysql
import sys
import math
import tempfile
import os.path
import subprocess
import difflib

db = _mysql.connect(host = "acbbdb1.picr.man.ac.uk", user = "nsmd", passwd = "nsmdP123", db = "nsmd")

if len(sys.argv) >= 2:
    workdir = sys.argv[1]
    workdirtemp = False
    sys.argv = sys.argv[1:]
else:
    workdir = tempfile.mkdtemp()
    workdirtemp = True

shortcodes = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "CSX": "c",
    "CSO": "d",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "MSE": "m",
    "PHE": "F",
    "PRO": "P",
    "SEP": "s",
    "SER": "S",
    "THR": "T",
    "TPO": "t",
    "TRP": "W",
    "TYR": "Y",
    "PTR": "y",
    "VAL": "V"
}

db = _mysql.connect(host = "acbbdb1.picr.man.ac.uk", user = "nsmd", passwd = "nsmdP123", db = "nsmd")

devnull = open("/dev/null", "w")

queryid = None
testpdb = None
accept_swaps = []
accept_misalign = False
accept_gap = True
verbose = True

for arg in sys.argv[1:]:
    if arg.startswith("accept="):
        subarg = arg[len("accept="):]
        if subarg == "misalign":
            accept_misalign = True
        elif subarg == "gap":
            accept_gap = True
        elif len(subarg) != 2:
            raise Exception("accept= parameters must either be two bases (e.g. accept=Ss) or accept=misalign or accept=gap")
        else:
            accept_swaps.append(subarg)
    elif arg.startswith("kinaseid="):
        queryid = int(arg[len("kinaseid="):])
    elif arg.startswith("testpdb="):
        testpdb = arg[len("testpdb="):]
    else:
        raise Exception("Bad argument %s" % arg)

query = "select gene, kinases.id, canonical_sequence, offset_canonical, length, pdbid, pdbchain, uniprotid from genes join kinases on gene_id = genes.id where "

if queryid is not None:
    query = "%s kinases.id = %d" % (query, queryid)
else:
    query = "%s pdbid is not null and pdbid != 'NOTFOUND'" % query

db.query(query)

if verbose:
    dbout = sys.stderr
else:
    dbout = open("/dev/null", "w")

for gene, kinaseid, trueseq, kinaseoffset, kinaselength, pdbid, pdbchain, uniprotid in db.store_result().fetch_row(maxrows=0):

    accept = True

    if testpdb is not None:
        testpdb = testpdb.split("/")
        pdbid = testpdb[0]
        pdbchain = testpdb[1]

    print >>dbout, "Check", gene, pdbid, pdbchain

    kinaseoffset = int(kinaseoffset)
    kinaselength = int(kinaselength)

    domainseq = trueseq[kinaseoffset : kinaseoffset + kinaselength]

    pdbfile = os.path.join(workdir, pdbid)

    if not os.path.exists(pdbfile):

        wget_cmd = ["/usr/bin/wget", "-O", pdbfile, "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=%s" % pdbid]
        try:
            subprocess.check_call(wget_cmd, stdout=devnull, stderr=devnull)
        except Exception as e:
            raise Exception("Fetch failed: command was " + " ".join(wget_cmd))

    pdbseq = ""

    startoff = -1
    endoff = -1

    other_ranges = []

    ignore_before = 0
    ignore_after = 0

    gaps = []

    pdb_up_offset = 0

    with open(pdbfile, "r") as f:
        for l in f:
            if l.startswith("SEQRES"):

                # We should have seen all the DBREFs and SEQADVs by now. Process them.

                for (start, end) in other_ranges:
                    if end < startoff:
                        ignore_before += ((end + 1) - start)
                    elif start > endoff:
                        ignore_after += ((end + 1) - start)
                    else:
                        print >>dbout, "\t", "DBREF range", start, end, "not handled"
                
                other_ranges = []

                chain = l[11]
                if chain != pdbchain:
                    continue
                bases = l[19:].split()
                for base in bases:
                    if ignore_before > 0:
                        ignore_before -= 1
                        continue
                    try:
                        if len(gaps) != 0 and (gaps[0][0] - 1) == len(pdbseq):
                            pdbseq = pdbseq + ("X" * (gaps[0][1] - gaps[0][0]))
                            gaps = gaps[1:]
                        pdbseq = pdbseq + shortcodes[base]
                    except Exception as e:
                        print >>dbout, "\t", "Reject base", base
                        pdbseq = pdbseq + "x"
                
            elif l.startswith("ATOM"):

                # We've seen all the SEQRES entries and therefore pdbseq is complete. Apply ignore_after.
                if ignore_after != 0:
                    pdbseq = pdbseq[:-ignore_after]
                    ignore_after = 0

                try:
                    offset = int(l[24:28]) + pdb_up_offset
                except:
                    print >>dbout, "\t", "Ignore atom with non-integer offset", l[24:28]
                chain = l[20:24].strip()
                base = shortcodes[l[17:20]]
                if chain != pdbchain:
                    continue

                if offset < startoff or offset > endoff:
                    continue
                elif (offset - 1) >= len(pdbseq):
                    print >>dbout, "\t", "Offset", offset, "out of range in ATOM record!"
                elif base != pdbseq[offset - 1]:
                    print >>dbout, "\t", "Disagreement in", pdbfile, "sequence has", pdbseq[offset - 1], "vs. atom line has", base, "at offset", offset

            elif l.startswith("DBREF"):

                bits = l.split()

                chain = bits[2]
                if chain != pdbchain:
                    continue

                db = bits[5]
                dbname = bits[6]
                if db != "UNP":
                    # These can indicate other blocks within the chain besides the uniprot mapping
                    #print >>dbout, "\t", "Other range", bits[-2], bits[-1]
                    other_ranges.append((int(bits[-2]), int(bits[-1])))
                    continue

                if dbname != uniprotid:
                    print >>dbout, "\t", "Warning! PDB file defines its offsets vs.", dbname, "not", uniprotid, "as expected"

                thisstart = int(bits[-2])
                localstart = int(bits[3])

                if startoff == -1:
                    startoff = thisstart
                    pdbseq = ("X" * (startoff - 1)) + pdbseq
                    pdb_up_offset = thisstart - localstart
                    if pdb_up_offset != 0:
                        print >>dbout, "Note: PDB vs. Uniprot coordinate offset", pdb_up_offset
                else:
                    gaps.append((endoff + 1, thisstart))
                    print >>dbout, "\t", "Gap (%d-%d]" % (endoff + 1, thisstart)
                    if localstart - thisstart != pdb_up_offset:
                        print >>dbout, "Warning! Multiple extents with different offsets vs. uniprot!"
       
                endoff = int(bits[-1])

            elif l.startswith("SEQADV"):

                bits = l.split()
                
                chain = bits[3]
                if chain != pdbchain:
                    continue

                offset = int(bits[4]) + pdb_up_offset
                if startoff == -1:
                    print >>dbout, "\t", "Ignored SEQADV because haven't seen DBREF"
                    continue
                if offset < startoff or offset > endoff:
                    #print >>dbout, "\t", "Other range (SEQADV)", offset, offset
                    other_ranges.append((offset, offset))

    pdbdomainseq = pdbseq[kinaseoffset : kinaseoffset + kinaselength]
                
    if startoff != -1:
        trueseq = trueseq[startoff - 1:]
        pdbseq = pdbseq[startoff - 1:]
        if kinaseoffset < (startoff - 1):
            domainseq = domainseq[startoff - (kinaseoffset + 1):]
            pdbdomainseq = pdbdomainseq[startoff - (kinaseoffset + 1):]
        trueseq = trueseq[:(endoff - startoff) + 1]
        pdbseq = pdbseq[:(endoff - startoff) + 1]
        if (kinaseoffset + kinaselength) > endoff:
            trimchars = (kinaseoffset + kinaselength) - endoff
            domainseq = domainseq[:-trimchars]
            pdbdomainseq = pdbdomainseq[:-trimchars]
    else:
        print >>dbout, "\t", "Canonical <-> PDB offsets not given?"

    if trueseq != pdbseq:
        #print >>dbout, "\t", "Whole sequence mismatch:"
        m = difflib.SequenceMatcher(None, trueseq, pdbseq, autojunk = False)

        #print >>dbout, "\t", trueseq
        #print >>dbout, "\t", pdbseq
        ops = m.get_opcodes()
        for (op, i1, i2, o1, o2) in ops:
            if op != 'equal':
                #print >>dbout, "\t", op, i1, i2, o1, o2, trueseq[i1:i2], pdbseq[o1:o2]
                pass
    if domainseq != pdbdomainseq:
        print >>dbout, "\t", "Kinase domain mismatch:"
        m = difflib.SequenceMatcher(None, domainseq, pdbdomainseq, autojunk = False)

        #print >>dbout, "\t", domainseq
        #print >>dbout, "\t", pdbdomainseq
        ops = m.get_opcodes()
        for (op, i1, i2, o1, o2) in ops:
            if op != 'equal':
                print >>dbout, "\t", op, i1, i2, o1, o2, domainseq[i1:i2], pdbdomainseq[o1:o2]
                if op == "replace" and i2 - i1 == o2 - o1:
                    charpairs = zip(domainseq[i1:i2], pdbdomainseq[o1:o2])
                    if all([(a + b) in accept_swaps or (b + a) in accept_swaps for (a, b) in charpairs]):
                        # Replacement consists only of allowed base replacements
                        continue
                    elif accept_gap and (domainseq[i1:i2] == ("X" * (i2 - i1)) or pdbdomainseq[o1:o2] == ("X" * (o2 - o1))):
                        # Replacement indicates a gap, nothing more
                        continue
                elif accept_misalign and (op == "insert" or op == "delete") and (i2 - i1 < 10 and o2 - o1 < 10) and (i1 == 0 or o1 == 0 or i2 == len(domainseq) or o2 == len(pdbdomainseq)):
                    # Short insertion or deletion at the start or finish of either sequence, indicating misalignment
                    continue
                # Otherwise the change is significant.
                accept = False

    if accept:
        print >>sys.stderr, "Accept", gene, pdbid, pdbchain
    else:
        print >>sys.stderr, "Reject", gene, pdbid, pdbchain

if workdirtemp:
    shutil.rmtree(workdir)
