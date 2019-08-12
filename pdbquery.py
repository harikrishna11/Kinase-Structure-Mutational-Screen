#!/usr/bin/python

import sys
import urllib2
import xml.etree.ElementTree as ET
import os.path
import subprocess
import socket
import random
import time

# Retries and timeouts are necessary because the PDB server has some amount of query rate-limiting, whose nature
# varies from month to month.

get_retries = 10
get_timeout = 3

def get(server, relurl, querydoc = None):

        try:
        
                for i in range(get_retries):
                        try:
                                req = urllib2.Request("http://%s%s" % (server, relurl), data=querydoc)
                                f = urllib2.urlopen(req, timeout = get_timeout)
                                return f.read()
                        except urllib2.URLError as e:
                                if not isinstance(e.reason, socket.timeout):
                                        raise e
                                print >>sys.stderr, "Timeout, retry..."
                        except socket.timeout as e:
                                print >>sys.stderr, "Timeout, retry..."

                raise Exception("Timed out %d times" % get_retries)

	except Exception as e:
		print >>sys.stderr, "Failed to retrieve", "http://%s%s" % (server, relurl), e
		print >>sys.stderr, "Query doc:"
		print >>sys.stderr, querydoc
		print >>sys.stderr, "\nResponse:"
		print >>sys.stderr, e.read()

def strip_suffix(pdbid):
        bits = pdbid.split(":")
        if len(bits) == 1 or len(bits) == 2:
                return bits[0]
        else:
                raise Exception("Unexpected PDB ID format: %s" % pdbid)

# Get PDB hits for the given genes, using query batching to speed the process up.
def getpdbinfo(names, verbose=False):

        gene_q = """
<?xml version="1.0" encoding="UTF-8"?>
<orgPdbQuery>
<queryType>org.pdb.query.simple.UpAccessionIdQuery</queryType>
<accessionIdList>%s</accessionIdList>
</orgPdbQuery>
""" % ",".join(names)

        codes = map(strip_suffix, get("www.rcsb.org", "/pdb/rest/search", gene_q).split())
        codes = filter(lambda x: x != 'null', codes)

        ret = dict()

        for name in names:
                ret[name] = dict()

        nonhuman_codes = 0
        parthuman_codes = 0
        unlabeled_codes = 0

        # For each PDB file returned:

        for code in codes:

                # Find resolution:
                if verbose:
                        print "Fetch", code
                descxml = get("www.rcsb.org", "/pdb/rest/describePDB?structureId=%s" % code)
                descroot = ET.fromstring(descxml)
                desc = descroot.find("PDB")
                try:
                        res = desc.attrib["resolution"]
                except KeyError:
                        # Acquired by method other than x-ray crystallography. Ignore.
                        continue
                try:
                        pdbname = desc.attrib["title"]
                except KeyError:
                        pdbname = "Untitled report %s" % code

                # Find ligands:
                ligxml = get("www.rcsb.org", "/pdb/rest/ligandInfo?structureId=%s" % code)
                ligroot = ET.fromstring(ligxml).find("ligandInfo")

                liglist = [lig.find("chemicalName").text for lig in ligroot.findall("ligand")]

                # Find relation to target genes:
                mapxml = get("www.rcsb.org", "/pdb/rest/das/pdb_uniprot_mapping/alignment?query=%s" % code)

                try:
                        maproot = ET.fromstring(mapxml)
                except:
                        print >>sys.stderr, "XML parse error for", code
                        continue

                for align in maproot.findall("{http://www.efamily.org.uk/xml/das/2004/06/17/dasalignment.xsd}alignment"):

                        uniqueobj = None

                        for alignobj in align.findall("{http://www.efamily.org.uk/xml/das/2004/06/17/dasalignment.xsd}alignObject"):

                                alignname = alignobj.attrib["dbAccessionId"]
                                if alignname not in names:
                                        continue
                                if uniqueobj is not None:
                                        print >>sys.stderr, "Object", code, "maps to at least", uniqueobj, "and", alignname
                                uniqueobj = alignname

                        if uniqueobj is None:
                                # This likely refers to a chain that is not related to any of the queried proteins
                                continue

                        unique_chain = None
                        regions = []
                        unique_offset = None

                        # Look for coverage records indicating that this PDB file covers a particular kinase domain.

                        for block in align.findall("{http://www.efamily.org.uk/xml/das/2004/06/17/dasalignment.xsd}block"):

                                # Expect to find two segment identifiers: one citing pdbid.pdbchain and the other citing name [a uniprotkb identifier]
                                segs = block.findall("{http://www.efamily.org.uk/xml/das/2004/06/17/dasalignment.xsd}segment")
                                if len(segs) != 2:
                                        print >>sys.stderr, "PDB entry", code, "has more than two segments in an alignment block"

                                pdbseg = None
                                uniprotseg = None

                                for seg in segs:
                                        if seg.attrib["intObjectId"] == uniqueobj:
                                                uniprotseg = seg
                                        elif seg.attrib["intObjectId"].startswith(code):
                                                pdbseg = seg

                                if pdbseg is None or uniprotseg is None:
                                        print >>sys.stderr, "Part of object", code, "could not be mapped onto uniprot entry", uniqueobj
                                        continue

                                try:

                                        chain = pdbseg.attrib["intObjectId"][len(code) + 1]
                                        region = (int(uniprotseg.attrib["start"]), int(uniprotseg.attrib["end"]))

                                        if unique_chain is not None and unique_chain != chain:
                                                raise Exception("PDB entry with multiple chains in one alignment")
                                        unique_chain = chain

                                        this_offset = int(uniprotseg.attrib["start"]) - int(pdbseg.attrib["start"])
                                        if unique_offset is None or unique_offset == this_offset:
                                                unique_offset = this_offset
                                        else:
                                                unique_offset = None

                                except ValueError:
                                        
                                        print >>sys.stderr, "Discard uniprot/pdb region with non-numeric coordinates"
                                        continue

                                regions.append(region)

                        if code not in ret[uniqueobj]:
                                ret[uniqueobj][code] = []
                        ret[uniqueobj][code].append({"resolution": float(res), "description": "%s chain %s (uniprot: %s)" % (pdbname, unique_chain, uniqueobj), "ligands": liglist, "chain": unique_chain, "regions": regions, "uniprot_pdb_offset": unique_offset})

        return ret

# Given a list of PDB entries retrieved using the query above, find the best one to model a given kinase domain.
# This checks that the PDB model in question is active or inactive also.

def getbestentry(name, results, seq, domain, workdir, uniprotid, find_active_model = True):

        domainoffset = seq.find(domain)
        if domainoffset == -1:
                print >>sys.stderr, name, "kinase domain not found in protein sequence"
                return None

        scored_results = []

        for code, entries in results.iteritems():
                
                for entry in entries:

                        # Filter out overlapping regions; some dubious PDB entries seem to have these.
                        # They may also occasionally have small backward alignments.

                        chain = entry["chain"]
                        regions = entry["regions"]

                        covered = []

                        for (start, stop) in regions:

                                if stop < start:
                                        continue

                                while len(covered) < stop:
                                        covered.append(False)

                                for i in range(start, stop + 1):
                                        covered[i - 1] = True # 1-based co-ordinates

                        # Attach a score to each entry: score 0 to 10 for domain coverage; score (3 - resolution) for res.
                        covered_bases = 0
                        for i in range(domainoffset, domainoffset + len(domain)):
                                if len(covered) > i and covered[i]:
                                        covered_bases += 1

                        scored_result = dict()
                        scored_result["coverage"] = (float(covered_bases) / len(domain))
                        scored_result["score"] = (10 * scored_result["coverage"]) + (3 - entry["resolution"])
                        scored_result["entry"] = entry
                        scored_result["chain"] = chain
                        scored_result["code"] = code
                        scored_result["uniprot_pdb_offset"] = entry["uniprot_pdb_offset"]
                        scored_results.append(scored_result)

        scored_results = sorted(scored_results, key = lambda x: x["score"], reverse=True)

        # Find an entry that we can confirm active
        while len(scored_results) > 0:
                
                pdb = scored_results[0]["code"]
                chain = scored_results[0]["chain"]

                pdbfname = os.path.join(workdir, pdb)
                if not os.path.exists(pdbfname):
                        print >>sys.stderr, "Fetch PDB file", pdb
                        pdbdata = get("www.rcsb.org", "/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=%s" % pdb)
                        with open(pdbfname, "w") as f:
                                f.write(pdbdata)

                with open(pdbfname, "r") as f:
                        lines = [l.strip() for l in f]
                        dbrefs = [l for l in lines if l.startswith("DBREF")]
                        remarks = [l for l in lines if l.startswith("REMARK   3 ")]

                # While we're looking at the PDB file, retrieve R-value stats that aren't available via the web interface.

                rvalue = None
                rfree = None
                for r in remarks:
                        bits = r[len("REMARK   3 "):].split()
                        if bits[:5] == ["R", "VALUE", "(WORKING", "SET)", ":"]:
                                rvalue = float(bits[5])
                        elif bits[:4] == ["FREE", "R", "VALUE", ":"]:
                                try:
                                        rfree = float(bits[4])
                                except ValueError as e:
                                        print >>sys.stderr, "Skip R-Free line", bits

                scored_results[0]["entry"]["rvalue"] = rvalue
                scored_results[0]["entry"]["rfree"] = rfree

                reject_mismatch = False

                # Check whether the PDB file covers the desired kinase domain.

                for dbref in dbrefs:
                        tag, pdbcode, pdbchain, pdbstart, pdbend, otherdb, otherdbid, otherdbdesc, otherstart, otherend = dbref.split()
                        if pdbchain != chain:
                                continue

                        if otherdb != "UNP":
                                continue

                        if otherdbid != uniprotid:

                                uniprotxml = get("www.uniprot.org", "/uniprot/%s.xml" % uniprotid)
                                xml = ET.fromstring(uniprotxml)
                                alternate_names = [a.text for a in xml.find("{http://uniprot.org/uniprot}entry").findall("{http://uniprot.org/uniprot}accession")]
                                if otherdbid in alternate_names:
                                        print >>sys.stderr, "Note: PDB references %s, which is an alias for the expected ID %s" % (otherdbid, uniprotid)
                                else:
                                        print >>sys.stderr, "Drop model", pdb, "that references", otherdbid, "not", uniprotid, "as expected"
                                        reject_mismatch = True
                                        break

                        given_offset = int(otherstart) - int(pdbstart)
                        actual_offset = scored_results[0]["uniprot_pdb_offset"]
                        if actual_offset is None:
                                print >>sys.stderr, "Drop model", pdb, "without simple linear mapping onto Uniprot", uniprotid
                                reject_mismatch = True
                                break
                        elif actual_offset != given_offset:
                                print >>sys.stderr, "Drop model", pdb, "whose stated uniprot-pdb offset", given_offset, "disagrees with the alignment offset", actual_offset, "(probably an isoform problem)"
                                reject_mismatch = True
                                break
             
                if reject_mismatch:
                        scored_results = scored_results[1:]
                        continue

                # Use check-active-dfg.py to determine whether the model is active/inactive as needed.

                p = subprocess.Popen(["./check-active-dfg.py", pdbfname, chain], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                (outs, errs) = p.communicate()
                bits = errs.split()
                if len(bits) != 0:
                        lastword = bits[-1]
                else:
                        lastword = None
                is_active = (lastword == "Active")
                is_inactive = (lastword == "Inactive")
                ret = p.returncode
                if ret != 0:
                        print >>sys.stderr, "Unexpected call result", ret, "stdout:", outs, "stderr:", errs
                        scored_results = scored_results[1:]                        
                elif (is_active and find_active_model) or (is_inactive and not find_active_model):
                        # Confirmed active/inactive, as required
                        break
                else:
                        print >>sys.stderr, "Drop", "inactive" if find_active_model else "active", "model", pdb, chain
                        scored_results = scored_results[1:]

        if len(scored_results) == 0:
                return None

        best_scored_result = scored_results[0]
        if best_scored_result["score"] < 5:
                print >>sys.stderr, "Drop entry with poor score", best_scored_result["code"], best_scored_result["chain"]
                return None

        ret = best_scored_result["entry"]
        ret["id"] = best_scored_result["code"]
        ret["best_chain"] = best_scored_result["chain"]
        ret["coverage"] = best_scored_result["coverage"]

        return ret
