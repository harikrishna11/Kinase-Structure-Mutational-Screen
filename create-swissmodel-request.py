#!/usr/bin/python

import sys
import _mysql
import json
import argparse

parser = argparse.ArgumentParser(description='Create a Swiss-Model request in JSON format')
parser.add_argument('output_file', metavar='outfile', help="Output JSON file")
parser.add_argument('--all', dest='all_files', action='store_true', default=False, help="Report all viable models")
parser.add_argument('--limit', metavar='N', dest='limit', type=int, default=-1, help="Report at most N models")
parser.add_argument('--match', metavar='PATTERN', dest='pattern', default="", help="Match mutation against SQL PATTERN")

args = vars(parser.parse_args())

if args["all_files"] == True and args["limit"] != -1:
    print >>sys.stderr, "--all and --limit are mutually exclusive"
    sys.exit(1)
elif args["all_files"] == False and args["limit"] == -1:
    print >>sys.stderr, "Must specify one of --all and --limit"
    sys.exit(1)

if args["pattern"] != "":
    pattern_string = "and sourceid like '%%%s%%'" % args["pattern"]
else:
    pattern_string = ""

db = _mysql.connect(host = "acbbdb1.picr.man.ac.uk", user = "nsmd", passwd = "nsmdP123", db = "nsmd")

query = "select sequences.id, sequence, pdbid, pdbchain from sequences inner join kinases on kinaseid = kinases.id where approved = TRUE and model_status = 0 %s order by kinaseid" % pattern_string

if args["limit"] != -1:
    query = query + " limit %d" % args["limit"]

db.query(query)
rows = db.store_result().fetch_row(maxrows=0)

with open(args["output_file"], "w") as f:
    
    outobj = [{"id": row[0], "sequence": row[1], "pdbid": row[2], "pdbchain": row[3]} for row in rows]
    json.dump(outobj, f)

update = "update sequences set model_status = 1 where id in (%s)" % ",".join([str(row[0]) for row in rows])
db.query(update)
db.store_result()

