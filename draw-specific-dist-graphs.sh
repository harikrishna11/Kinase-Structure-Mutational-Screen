#!/bin/bash

wildpath=$2

MYDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

for f in `ls $1/*_to_*.xvg`; do
$MYDIR/plot-mavg-lines.py "${f/.xvg/} inter-atomic distance vs. time" "Time (ns)" "Distance (nm)" $f $f $wildpath/$f WT $f.pdf
done
