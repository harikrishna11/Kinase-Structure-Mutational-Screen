#!/bin/bash

MYDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

mkdir -p $1/joint_analysis
python $MYDIR/gather-specific-graphs.py $1 $1/joint_analysis
cd $1/joint_analysis
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=merged.pdf `ls *.pdf | sort`
