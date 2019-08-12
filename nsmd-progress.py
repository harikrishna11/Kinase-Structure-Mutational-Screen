#!/usr/bin/env python

import os
import time
import sys
import subprocess

rootdir = "/data/snc/nstephenson/nsmd/results"
scriptdir = "/data/snc/nstephenson/nsmd/scripts"

inprogress_reports = []
analysis_inprogress_reports = []
done_reports = []
unknown_reports = []

phases = ["main_sim_400ns", "NPT", "NVT", "EM"]
phase_mdps = ["sim", "npt", "nvt", "minim"]

def read_nsteps(phasename):
    mdp = os.path.join(scriptdir, "%s.mdp" % phasename)
    with open(mdp, "r") as f:
        for l in f:
            bits = l.split()
            if len(bits) >= 3 and bits[0] == "nsteps" and bits[1] == "=":
                return int(bits[2])
    raise Exception("No nsteps in %s" % mdp)

def line_is_progressheader(l):
    bits = tuple(l.strip().split())
    return bits == ("Step", "Time", "Lambda")

def read_nsteps_progress(l):
    bits = tuple(l.strip().split())
    return int(bits[0])

phase_steps = map(read_nsteps, phase_mdps)

def check_jobfile(jobfile):

    with open(jobfile, "r") as f:
        oldjobid = f.read().strip()
        with open("/dev/null", "w") as nullout:
            checkret = subprocess.Popen(["/apps/modules/pkg/clusterware/torque/5.1.0-1/gcc-4.4.7/bin/qstat", oldjobid], stdout=nullout, stderr=subprocess.STDOUT).wait()
        return checkret == 0

for gene in os.listdir(rootdir):

    if gene == "ABL1_MODELS":
        continue
    
    for runspec in os.listdir(os.path.join(rootdir, gene)):

        if runspec == "joint_analysis":
            continue

        runpath = os.path.join(rootdir, gene, runspec)
        if not os.path.isdir(runpath):
            continue

        jobfile = os.path.join(runpath, "jobid.txt")
        analysis_jobfile = os.path.join(runpath, "analysis_jobid.txt")

        running = False
        analysing = False
        if os.path.exists(jobfile):
            running = check_jobfile(jobfile)
        if os.path.exists(analysis_jobfile):
            analysing = check_jobfile(analysis_jobfile)
                    
        final_donefile = os.path.join(runpath, "main_sim_400ns_done.flag")
        analysis_donefile = os.path.join(runpath, "400ns_specific_graphs.done")
        if os.path.exists(final_donefile):
            analysis_done = os.path.exists(analysis_donefile)
            if (not analysis_done) and analysing:
                analysis_inprogress_reports.append((gene, runspec))
            else:
                done_reports.append((gene, runspec, os.stat(final_donefile), running, analysis_done))
        else:
            report = None
            for (p, steps) in zip(phases, phase_steps):
                logfile = os.path.join(runpath, "%s.log" % p)
                donefile = os.path.join(runpath, "%s_done.flag" % p)
                stepsdone = None
                if os.path.exists(logfile) and not os.path.exists(donefile):
                    with open(logfile, "r") as f:
                        f.seek(-10000, 2)
                        ls = f.readlines()
                        for i in range(len(ls) - 1, -1, -1):
                            if line_is_progressheader(ls[i]):
                                stepsdone = read_nsteps_progress(ls[i+1])
                                break
                if stepsdone is not None:
                    report = (gene, runspec, p, float(stepsdone) / steps, running)
                    break
            if report is not None:
                inprogress_reports.append(report)
            else:
                unknown_reports.append((gene, runspec, running))

def pp_running(running):
    return "RUNNING" if running else "STOPPED"

inprogress_reports = sorted(inprogress_reports)
unknown_reports = sorted(unknown_reports)
done_reports = sorted(done_reports, key = lambda x : x[2].st_mtime, reverse = True)

categories = False

outlines = []

def header(ln):
    if categories:
        outlines.append(ln)

if len(inprogress_reports) != 0:
    header("In progress:")
for (gene, runspec, phasename, doneprop, running) in inprogress_reports:
    outlines.append("  %s / %s: %s (%.1f%% done) [%s]" % (gene, runspec, phasename, doneprop * 100, pp_running(running)))
if len(analysis_inprogress_reports) != 0:
    header("Analysis in progress:")
for (gene, runspec) in analysis_inprogress_reports:
    outlines.append("  %s / %s: analysis in progress" % (gene, runspec))
if len(unknown_reports) != 0:
    header("Unknown state:")
for (gene, runspec, running) in unknown_reports:
    outlines.append("  %s / %s: unknown state [%s]" % (gene, runspec, pp_running(running)))
if len(done_reports) != 0:
    header("Done:")
for (gene, runspec, stat, running, analysis_done) in done_reports:
    outlines.append("  %s / %s: done (completed %s) [%s]" % (gene, runspec, time.ctime(stat.st_mtime), "Analysis complete" if analysis_done else "Analysis not started"))

if not categories:
    outlines = sorted(outlines)

for l in outlines:
    print l
        
