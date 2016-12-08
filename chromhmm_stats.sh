#!/usr/bin/env bash

# Script to analyze chromhmm states composition of given file
# author: oleg.shpynov@jetbrains.com

# Example:
# > chromhmm_states.sh k27ac.diff cd14_chromhmm.bed
# 12_EnhBiv	200
# 13_ReprPC	400
# 15_Quies	800
# 4_Tx	800

# Check tool.
which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }

if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <BED_FILE> <CHROMHMM_MARKUP_FILE>"
    exit 1
fi

FILE=$1
CHROMHMM=$2

# Compute intersection length by chromhmm state
bedtools intersect -a $FILE -b $CHROMHMM -wa -wb |\
 awk -v COLS=$(head -1 $FILE | awk '{ print NF }') '{ tot[$(4+COLS)]+=(($(3+COLS)<$3?$(3+COLS):$3)-($(2+COLS)>$2?$(2+COLS):$2)) } END { for (i in tot) print i"\t"tot[i] }' |\
sort