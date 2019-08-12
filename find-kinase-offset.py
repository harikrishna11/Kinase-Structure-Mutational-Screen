
import sys
import _mysql
import difflib

db = _mysql.connect(host = "acbbdb1.picr.man.ac.uk", user = "nsmd", passwd = "nsmdP123", db = "nsmd")

def getuniprotmatch(gene):

    db.query("select id, matched_sequence, kinase_offset from genes where gene = '%s'" % gene)
    rows = db.store_result().fetch_row(maxrows=0)
    if len(rows) == 1:
        return rows[0]
    elif len(rows) > 1:
        raise Exception("Multiple uniprot IDs for %s" % gene)
    else:
        return (None, None, None)

db.query("select gene, given_sequence, matched_sequence, kinases.id, offset_given, length from genes inner join kinases on genes.id = kinases.gene_id where genes.matched_sequence is not null and offset_matched is null")
for gene, given_seq, matched_seq, record_id, kinase_offset, kinase_length in db.store_result().fetch_row(maxrows=0):
              
    kinase_offset = int(kinase_offset)
    kinase_length = int(kinase_length)

    kinase_end = kinase_offset + kinase_length
    new_kinase_offset = kinase_offset

    m = difflib.SequenceMatcher(None, given_seq, matched_seq, autojunk = False)

    # Kinase is at kinase_offset in left-hand sequence.
    # Find out where it ends up, and note any length-preserving substitutions that affect it.

    ops = m.get_opcodes()
    for tag, i1, i2, j1, j2 in ops:

        if tag == "equal":
            # No effect
            continue

        lenchg = (j2 - j1) - (i2 - i1)

        if i2 <= kinase_offset:
            # Alteration before kinase. May move the kinase.
            new_kinase_offset += lenchg

        elif i1 >= kinase_end:
            # Alteration after kinase. No effect.
            continue

        else:
            # Alteration partially or wholly ovelaps the kinase.
            # Reject for manual inspection if it changes length:
            if lenchg != 0:
                print >>sys.stderr, "Skip", gene, "(Length-changing alteration overlaps kinase domain)"
                new_kinase_offset = None
                break

            # Otherwise note that the kinase is altered.
            if i1 < kinase_offset:
                chg = kinase_offset - i1
                i1 += chg
                j1 += chg
            if i2 > kinase_end:
                chg = i2 - kinase_end
                i2 -= chg
                j2 -= chg

            print >>sys.stderr, "Warning: kinase domain", record_id, "altered:", given_seq[i1:i2], "->", matched_seq[j1:j2], "at offset", i1

    if new_kinase_offset is None:
        # Manual intervention required
        continue
    if new_kinase_offset != kinase_offset:
        print >>sys.stderr, "Warning: kinase domain", record_id, "moved on mapping to reference sequence (offset %d)" % (new_kinase_offset - kinase_offset)

    db.query("update kinases set offset_matched = %d where id = %s" % (new_kinase_offset, record_id))
    db.store_result()

