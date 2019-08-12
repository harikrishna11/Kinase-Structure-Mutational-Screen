#!/bin/bash

wildpath=$2

MYDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

for f in `ls $1/*.xvg | grep -v _to_`; do
$MYDIR/plot-mavg-lines.py "RMSD vs. time for ${f/.xvg/}" "Time (ns)" "RMSD" $f $f $wildpath/$f WT $f.pdf
done
