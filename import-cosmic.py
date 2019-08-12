#!/usr/bin/python

import _mysql
import sys
import mutate

if len(sys.argv) < 2:
    print >>sys.stderr, "Usage: import-cosmic.py datafile.tsv"

db = _mysql.connect(host = "acbbdb1", user = "nsmd", passwd = "nsmdP123", db = "nsmd")

# Fetch list of kinase domains for each gene

db.query("select gene, offset_canonical, length from genes join kinases on gene_id = genes.id where pdbid is not null and pdbid != 'NOTFOUND'")
rows = db.store_result().fetch_row(maxrows=0)

gene_ranges = dict()

def aa_in_scope(gene, offset):
    try:
        ranges = gene_ranges[gene]
    except KeyError:
        return None

    for (roffset, rlength) in ranges:
        if offset >= roffset and offset < (roffset + rlength):
            return (roffset, rlength)
    return None

for (gene, offset, length) in rows:
    
    if gene not in gene_ranges:
        gene_ranges[gene] = []
    gene_ranges[gene].append((int(offset), int(length)))

# Fetch canonical sequence for each gene of interest

db.query("select gene, canonical_sequence from genes")
rows = db.store_result().fetch_row(maxrows=0)

seqs = dict()

for (gene, sequence) in rows:
    seqs[gene] = sequence

# Organise the dataset by cosmic identifier.

recs_by_cosm = dict()

with open(sys.argv[1], "r") as f:

    for rec in f:

        fields = [f.strip() for f in rec.split("\t")]
        cosm = fields[0]
        if cosm not in recs_by_cosm:
            recs_by_cosm[cosm] = []
        recs_by_cosm[cosm].append(fields)

# For genes with multiple transcripts mentioned in the dataset, find which one aligns best with our canonical sequence.

alignments_by_gene = dict()

for recs in recs_by_cosm.itervalues():
    
    for rec in recs:

        thisgene = rec[9]
        thistranscript = rec[8]
        offset = int(rec[6])
        aa_change = rec[7]

        try:
            base_seq = seqs[thisgene]
            gene = thisgene
        except KeyError:
            continue

        (oldaa, newaa) = aa_change.split("/")

        if thisgene not in alignments_by_gene:
            alignments_by_gene[thisgene] = dict()
        if thistranscript not in alignments_by_gene[thisgene]:
            alignments_by_gene[thisgene][thistranscript] = [0, 0]
        alignments_by_gene[thisgene][thistranscript][1] += 1
        if len(base_seq) > (offset - 1) and base_seq[offset - 1] == oldaa:
            alignments_by_gene[thisgene][thistranscript][0] += 1

best_transcript_by_gene = dict()

for gene, transcripts in alignments_by_gene.iteritems():

    transcripts = sorted(list(transcripts.items()), key = lambda x : x[1], reverse = True)
    if len(transcripts) > 1 and transcripts[0][1][0] == transcripts[1][1][0]:
        print "Gene", gene, "transcripts", transcripts[0][0], transcripts[1][0], "tied for first place (%d/%d each)" % (transcripts[0][1][0], transcripts[0][1][1])
    if transcripts[0][1][0] != transcripts[0][1][1]:
        print "Best transcript for", gene, transcripts[0][0], "only hit %d/%d times" % (transcripts[0][1][0], transcripts[0][1][1])
    else:
        print "Best transcript for", gene, "hit perfectly (%d)" % transcripts[0][1][0]

    best_transcript_by_gene[gene] = transcripts[0][0]

nomatch_by_gene = dict()

stats_by_gene = dict()

for recs in recs_by_cosm.itervalues():

    matching_rec = None
    gene = None
    failures = []
    in_scope = None
    mutated_seq = None
    for rec in recs:
        
        thisgene = rec[9]
        thistranscript = rec[8]
        offset = int(rec[6])
        aa_change = rec[7]

        try:
            base_seq = seqs[thisgene]
            gene = thisgene
        except KeyError:
            continue
        
        if thistranscript != best_transcript_by_gene[gene]:
            continue

        (oldaa, newaa) = aa_change.split("/")
        in_scope = aa_in_scope(thisgene, offset - 1)
        try:
            mutated_seq = mutate.applymutation(base_seq, oldaa, newaa, offset)
            if matching_rec is not None:
                raise Exception("Two records with same transcript? %s, %s" % (rec, matching_rec))
            matching_rec = rec
        except Exception as e:
            failures.append(e)

    if gene not in stats_by_gene:
        stats_by_gene[gene] = {"nomatch": 0, "nomatch_kinase": 0, "match": 0, "match_kinase": 0, "match_kinase_new": 0}

    if gene is None:
        pass
    elif matching_rec is None:
        stats_by_gene[gene]["nomatch"] += 1
        if in_scope is not None:
            stats_by_gene[gene]["nomatch_kinase"] += 1
    else:
        stats_by_gene[gene]["match"] += 1
        if in_scope is not None:
            stats_by_gene[gene]["match_kinase"] += 1
            kinase_offset, kinase_length = in_scope
            db.query("select id from sequences where sequence = '%s'" % mutated_seq[kinase_offset : kinase_offset + kinase_length])
            rows = db.store_result().fetch_row(maxrows=0)
            if len(rows) == 0:
                stats_by_gene[gene]["match_kinase_new"] += 1
            
stats = ["match", "match_kinase", "match_kinase_new", "nomatch", "nomatch_kinase"]

for stat in stats:
    print "Total", stat, sum([stats[stat] for stats in stats_by_gene.itervalues()])

for stat in stats:
    print "Top 10", stat
    tops = sorted(list(stats_by_gene.iteritems()), key = lambda x: x[1][stat], reverse = True)
    for gene, statblock in tops[:10]:
        print gene, statblock[stat]

