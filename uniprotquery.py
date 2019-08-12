
import urllib2
import xml.etree.ElementTree as ET
import sys
import difflib

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

def apply_modification(modseq, modtag):

    original = modtag.find("{http://uniprot.org/uniprot}original")
    variation = modtag.find("{http://uniprot.org/uniprot}variation")
    loctag = modtag.find("{http://uniprot.org/uniprot}location")
    postag = None
    begintag = None
    endtag = None
    if loctag is not None:
        postag = loctag.find("{http://uniprot.org/uniprot}position")
        begintag = loctag.find("{http://uniprot.org/uniprot}begin")
        endtag = loctag.find("{http://uniprot.org/uniprot}end")

    begin = None
    end = None

    # 1-based, inclusive coordinates converted to string coordinates.

    if postag is not None:
        begin = int(postag.attrib["position"]) - 1
        end = begin + 1
    elif begintag is not None and endtag is not None:
        begin = int(begintag.attrib["position"]) - 1
        end = int(endtag.attrib["position"])

    handled = False

    if begin is not None and end is not None:

        if original is not None and variation is not None:

            origbase = original.text
            varbase = variation.text

            if len(origbase) != (end - begin):
                print >>sys.stderr, "Mutation spec", ET.tostring(modtag), "lengths inconsistent"
                return

            # Substitution?

            for i in range(begin, end):

                if modseq[i] != origbase[i - begin]:
                    print >>sys.stderr, "Mutation spec", ET.tostring(modtag), "clashes at position", i
                    return
                modseq[i] = ""

            if len(varbase) == len(origbase):

                # If not length-modifying, overwrite each base:

                for i in range(begin, end):
                    modseq[i] = varbase[i - begin]

            else:
                
                # Otherwise represent like a large insertion:
                modseq[begin] = varbase

            handled = True

        else:

            # Deletion specifier?
            handled = True

            for i in range(begin, end):
                modseq[i] = ""

    if not handled:
        print >>sys.stderr, "Modification", ET.tostring(modtag), "not applied"

def get_opcode_string(matcher, s1, s2):

    instructions = []

    for tag, i1, i2, j1, j2 in matcher.get_opcodes():

        if tag == "equal":
            instructions.append("%d equal bases" % (i2 - i1))
        elif tag == "delete":
            instructions.append("delete %s" % s1[i1:i2])
        elif tag == "insert":
            instructions.append("insert %s" % s2[j1:j2])
        else:
            instructions.append("replace %s with %s" % (s1[i1:i2], s2[j1:j2]))

    return ", ".join(instructions)

def getisoforms(result, baseseq, resname):

    altprods = None

    for comment in result.findall("{http://uniprot.org/uniprot}comment"):
        if comment.attrib["type"] == "alternative products":
            altprods = comment

    if altprods is None:
        return {resname: baseseq}

    seqs = dict()
    mods = dict()
    for feature in result.findall("{http://uniprot.org/uniprot}feature"):
        if "id" in feature.attrib:
            mods[feature.attrib["id"]] = feature

    for isoform in altprods.findall("{http://uniprot.org/uniprot}isoform"):

        newseqname = isoform.find("{http://uniprot.org/uniprot}id").text

        newseqtag = isoform.find("{http://uniprot.org/uniprot}sequence")
        if newseqtag is None:
            print >>sys.stderr, "isoform without sequence processing %s?" % gene
            continue

        if newseqtag.attrib["type"] == "described":

            modseq = []
            for i in range(len(baseseq)):
                modseq.append(baseseq[i])

            for applymod in newseqtag.attrib["ref"].split():

                try:
                    modtag = mods[applymod]
                except KeyError:
                    print >>sys.stderr, "No such mod", applymod, "processing isoforms for", gene
                    continue

                apply_modification(modseq, modtag)

            seqs[newseqname] = ("".join(modseq))

        elif newseqtag.attrib["type"] == "displayed":

            seqs[newseqname] = baseseq

        else:

            print >>sys.stderr, "Unknown isoform tag", ET.tostring(isoform)

    return seqs
    
def getbestmatch(gene, seq):

    # Hack: TTN takes an age to analyse
    if gene == "TTN":
        return {"uniprotid": "Q8WZ42", "isoform": "Q8WZ42-1", "score": 1.0, "diffops": "%d equal bases" % len(seq), "matched_sequence": seq}

    results = get("www.uniprot.org", "/uniprot/?query=gene%%3a%s+AND+organism%%3a\"Homo+sapiens\"&format=xml" % gene)
    if len(results) == 0:
        return None

    try:
        resultroot = ET.fromstring(results)
    except:
        print >>sys.stderr, "XML parse error for uniprot query", gene
        return None

    bestresult = None
    bestmatch = 0.0

    for result in resultroot.findall("{http://uniprot.org/uniprot}entry"):

        resname = result.find("{http://uniprot.org/uniprot}accession").text

        baseseq = result.find("{http://uniprot.org/uniprot}sequence")
        if baseseq is None:
            continue
        baseseq = "".join(baseseq.text.split())

        seqs = getisoforms(result, baseseq, resname)

        for seqname, rseq in seqs.iteritems():

            s = difflib.SequenceMatcher(None, seq, rseq, autojunk = False)
            match = s.ratio()
            print >>sys.stderr, seqname, match
            if match > bestmatch:
                bestresult = {"uniprotid": resname, "isoform": seqname, "score": match, "diffops": get_opcode_string(s, seq, rseq), "matched_sequence": rseq}
                bestmatch = match

    if bestresult is not None and bestmatch < 0.9:
        print >>sys.stderr, "Drop best match (%g) for %s" % (bestmatch, gene)
        bestresult = None
    elif bestmatch < 1.0:
        print >>sys.stderr, "Using inexact match (match %g) for %s" % (bestmatch, gene)

    return bestresult

