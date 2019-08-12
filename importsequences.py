#!/usr/bin/python

import sys
import _mysql
import difflib

def getgeneid(db, gene):

        existing_gene = db.query("select id from genes where gene = '%s'" % gene)
        rows = db.store_result().fetch_row(maxrows=0)
        if len(rows) > 1:
                raise Exception("Duplicates for %s" % gene)
        elif len(rows) == 1:
                return rows[0][0]
        else:
                return None

def get_kinase_offset(seq, kinase):
        if kinase == "NULL":
                return None
        offset = seq.find(kinase)
        if offset == -1:
                qual = ""
                chunks = 8
                for i in range(chunks):
                        chunk = len(kinase) / chunks
                        off = chunk * i
                        subseq = kinase[off:off+chunk]
                        found = seq.find(subseq)
                        if found != 0 and (found - off) >= 0:
                                try_kinase = seq[(found - off) : (found - off) + len(kinase)]
                                if len(try_kinase) != len(kinase):
                                        continue
                                m = difflib.SequenceMatcher(None, kinase, try_kinase, autojunk = False)
                                if m.ratio() > 0.9:
                                        print >>sys.stderr, "Using inexact kinase match"
                                        return (found - off)

                raise Exception("Kinase not found in given sequence%s" % qual)
        return offset

def import_gene(db, gene, sequence, kinases):

                geneid = getgeneid(db, gene)
                
                if geneid is not None:
                        print >>sys.stderr, "Skip", gene, "already in DB"
                        return
                
                try:
                        kinase_offsets = [get_kinase_offset(sequence, x) for x in kinases]
                except Exception as e:
                        print >>sys.stderr, "Skip", gene, e
                        return
                        
                db.query("insert into genes (gene, given_sequence) values ('%s', '%s')" % (gene, sequence))
                db.store_result()

                geneid = int(getgeneid(db, gene))

                for kinase, kinase_offset in zip(kinases, kinase_offsets):
          
                        if kinase_offset is None:
                                continue
                        db.query("insert into kinases (gene_id, offset_given, length) values (%d, %d, %d)" % (geneid, kinase_offset, len(kinase)))

