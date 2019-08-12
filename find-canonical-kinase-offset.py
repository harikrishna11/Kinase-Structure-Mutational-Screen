
import urllib2
import xml.etree.ElementTree as ET
import sys
import uniprotquery
import _mysql

db = _mysql.connect(host = "acbbdb1.picr.man.ac.uk", user = "nsmd", passwd = "nsmdP123", db = "nsmd")

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

def map_to_canonical_sequence(uniprot_id, isoform_id, kinase_offset, kinase, gene):

    results = get("www.uniprot.org", "/uniprot/%s.xml" % uniprot_id)
    result = ET.fromstring(results).find("{http://uniprot.org/uniprot}entry")

    altprods = None

    print uniprot_id, isoform_id

    for comment in result.findall("{http://uniprot.org/uniprot}comment"):
        if comment.attrib["type"] == "alternative products":
            altprods = comment

    # Find the modifications that will be used to describe the isoform:
    mods = dict()
    for feature in result.findall("{http://uniprot.org/uniprot}feature"):
        if "id" in feature.attrib:
            mods[feature.attrib["id"]] = feature

    
    # Find the requested isoform
    isoform = None

    for thisisoform in altprods.findall("{http://uniprot.org/uniprot}isoform"):
        if thisisoform.find("{http://uniprot.org/uniprot}id").text == isoform_id:
            isoform = thisisoform

    baseseq = result.find("{http://uniprot.org/uniprot}sequence")
    baseseq = "".join(baseseq.text.split())

    modseq = []
    for i in range(len(baseseq)):
        modseq.append(baseseq[i])

    seqtag = isoform.find("{http://uniprot.org/uniprot}sequence")
    if seqtag.attrib["type"] == "displayed":
        # This is the canonical sequence!
        return (kinase_offset, baseseq)

    # Otherwise, map from a noncanonical sequence to the canonical one.

    for applymod in isoform.find("{http://uniprot.org/uniprot}sequence").attrib["ref"].split():

        modtag = mods[applymod]
        uniprotquery.apply_modification(modseq, modtag)

    # Find the old offset (i.e. offset in the canonical string) that gives rise to kinase_offset in the isoform:
    kinase_end = kinase_offset + len(kinase)
    new_kinase_offset = None
    c = 0

    for i in range(len(modseq)):

        if new_kinase_offset is None:

            if c > kinase_offset:
                print >>sys.stderr, "Splice overlapping start of kinase for", gene, isoform_id
                return (None, None)
            elif c == kinase_offset:
                # Kinase begins here before splicing
                new_kinase_offset = i
            # Else we're not up the kinase yet

        else:

            if c == kinase_end:
                break
            if len(modseq[i]) != 1:
                print >>sys.stderr, "Length-changing alteration in kinase for", gene, isoform_id
                return (None, None)
                
        c += len(modseq[i])

    if new_kinase_offset is None:
        print >>sys.stderr, "Never found start of kinase for", gene, isoform_id
        return (None, None)

    return new_kinase_offset, baseseq

def getuniprotmatch(gene):

    db.query("select id, uniprotid, isoform, matched_sequence, kinase_offset, kinase_offset_canonical from genes where gene = '%s'" % gene)
    rows = db.store_result().fetch_row(maxrows=0)
    if len(rows) == 1:
        return rows[0]
    elif len(rows) > 1:
        raise Exception("Multiple uniprot IDs for %s" % gene)
    else:
        return (None, None, None, None, None, None)

db.query("select gene, genes.id, kinases.id, uniprotid, matched_isoform, matched_sequence, offset_matched, length from genes inner join kinases on genes.id = kinases.gene_id where offset_matched is not null and offset_canonical is null")
for gene, gene_id, record_id, uniprot_id, isoform, matched_seq, matched_offset, kinase_length in db.store_result().fetch_row(maxrows=0):

    matched_offset = int(matched_offset)
    kinase = matched_seq[matched_offset : matched_offset + int(kinase_length)]

    if isoform == uniprot_id:

        # Already mapped to the correct isoform
        new_kinase_offset = matched_offset
        kinase_moved = 0
        baseseq = matched_seq

    else:

        new_kinase_offset, baseseq = map_to_canonical_sequence(uniprot_id, isoform, matched_offset, kinase, gene)
        if new_kinase_offset is None:
            continue

    db.query("update kinases set offset_canonical = %d where id = %s" % (new_kinase_offset, record_id))
    db.store_result()

    db.query("update genes set canonical_sequence = '%s' where id = %s" % (baseseq, gene_id))
    db.store_result()

