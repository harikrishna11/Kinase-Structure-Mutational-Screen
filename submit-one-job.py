
import _mysql
import PBSPy.capi as pbs
import os.path
import sys

db = _mysql.connect(host = "acbbdb1.picr.man.ac.uk", user = "nsmd", passwd = "nsmdP123", db = "nsmd")
torque = pbs.Server()
torque.connect()

data_root = "/data/snc/nstephenson/nsmd/results"
qname = "bioinf.q"
cores = 8

query = "select id from sequences where model_status = 2"
if len(sys.argv) > 1:
    query = "%s and id = %d" % (query, int(sys.argv[1]))

# Find a candidate that isn't currently running:
db.query(query)
rows = db.store_result().fetch_row(maxrows=0)

runnable_ids = dict()
for row in rows:
    runnable_ids[int(row[0])] = set([1, 2, 3])
    
# Eliminate jobs that are currently running or queued against the cluster:
jobs = torque.selectjob()
for job in jobs:
    attrs = torque.statjob(job)[0].attribs
    for a in attrs:
        if a.name == "Job_Name" and a.value.startswith("nsmd"):
            bits = a.value.split(".")
            thisid = int(bits[1])
            thisrun = int(bits[2])
            if thisid not in runnable_ids:
                print >>sys.stderr, "Anomaly: job %s running or scheduled but not in runnable_ids" % a.value
                break
            if thisrun not in runnable_ids[thisid]:
                print >>sys.stderr, "Anomaly: job %s running or scheduled but not in runnable_ids[%d]" % (a.value, thisid)
                break
            runnable_ids[thisid].remove(thisrun)
            if len(runnable_ids[thisid]) == 0:
                del runnable_ids[thisid]

todel = []

# Eliminate jobs that have completion flags set:
for thisid, runs in runnable_ids.iteritems():
    for run in runs:
        path = os.path.join(data_root, str(thisid), "gmx-%d" % run, "main_sim_done.flag")
        if os.path.exists(path):
            todel.append((thisid, run))

for thisid, run in todel:
    runnable_ids[thisid].remove(run)
    if len(runnable_ids[thisid]) == 0:
        del runnable_ids[thisid]

runnable_ids = sorted(runnable_ids.iteritems(), reverse = True, key = lambda kv: len(kv[1]))

for runid, runs in runnable_ids:
    for run in runs:
        print "Will run", runid, run

        makefile = os.path.join(data_root, "Makefile")
        workdir = os.path.join(data_root, str(runid), "gmx-%d" % run)
        script = os.path.join(workdir, "run.torque")
        stdout = os.path.join(workdir, "nsmd.stdout.log")
        stderr = os.path.join(workdir, "nsmd.stderr.log")

        try:
            os.mkdir(workdir)
        except OSError as e:
            if e.errno == 17:
                print >>sys.stderr, "Warning", workdir, "already exists; looks like a resume?"
            else:
                raise e

        with open(script, "w") as f:
            f.write("/usr/bin/make -f %s -C %s" % (makefile, workdir))

        attribs = [
            pbs.Attr("Output_Path", stdout),
            pbs.Attr("Error_Path", stderr),
            pbs.Attr("Job_Name", "nsmd.%d.%d" % (runid, run)), 
            pbs.Attr("Resource_List", "1:ppn=%d" % cores, "nodes")
        ]

        newjob = torque.submit(attribs, script)
        print "Started as", newjob
        
        sys.exit(0)
