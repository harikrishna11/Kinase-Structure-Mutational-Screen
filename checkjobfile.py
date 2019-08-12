
import os.path
import subprocess
import sys

def checkjobfile(jobfile):

    if os.path.exists(jobfile):
        with open(jobfile, "r") as f:
            oldjobid = f.read().strip()
        with open("/dev/null", "w") as nullout:
            checkret = subprocess.Popen(["/apps/modules/pkg/clusterware/torque/5.1.0-1/gcc-4.4.7/bin/qstat", oldjobid], stdout=nullout, stderr=subprocess.STDOUT).wait()
        if checkret == 0:
            raise Exception("This run seems to already be in progress, as job %s" % oldjobid)
        elif checkret != 153:
            raise Exception("Unexpected return code %d from qstat %s" % (checkret, oldjobid))
        else:
            print >>sys.stderr, "This run seems to have existed as job %s in the past, but that job no longer exists; will restart..." % oldjobid
