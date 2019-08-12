#!/bin/bash

if [ $# -ne 4 ]; then
    echo "Usage: manual-rmsd-analysis.sh reference-directory trajectory-directory output-directory output-prefix"
    exit 1
fi

MYDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

mkdir -p $3
jobname=`printf "%s_vs_%s" \`basename $1\` \`basename $2\``
echo $MYDIR/run-rmsd-analyses.py $3 $4 gmx rms -s `readlink -f $1/NPT_CA.gro` -f `readlink -f $2/main_sim_400ns_nopbc_CA.xtc` -n `readlink -f $2/main_sim_400ns_ca_index.ndx` | qsub -o $3/stdout.log -e $3/stderr.log -N $jobname -m a -M Natalie.Stephenson@cruk.manchester.ac.uk -lwalltime=20:00:00:00


