#!/usr/bin/env bash

# Script to analyze chromhmm states composition of given file
# author: oleg.shpynov@jetbrains.com
#
# ChromHMM precomputed markup is available at:
# https://www.encodeproject.org/report/?type=Annotation&software_used.software.name=chromhmm

# Example:
# > chromhmm_states.sh k27ac.diff cd14_chromhmm.bed
# 12_EnhBiv	200
# 13_ReprPC	400
# 15_Quies	800
# 4_Tx	800

# Check tool.
which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }

if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <BED_FILE> <CHROMHMM_MARKUP_FILE.BED>"
    exit 1
fi
>&2 echo "chromhmm_stat $@"

FILE=$1
CHROMHMM=$2

if [ ! $CHROMHMM =~ *.bed ]; then
    which bigBedToBed &>/dev/null || { echo "bigBedToBed not found! Download: http://hgdownload.cse.ucsc.edu/admin/exe/"; exit 1; }
    NAME=${CHROMHMM%%.*}
    echo "Bed file required, found bigBed. Converting to $CHROMHMM.bed"
    bigBedToBed $FILE $CHROMHMM.bed
    CHROMHMM=$CHROMHMM.bed
fi

# Compute intersection length by chromhmm state
bedtools intersect -a $FILE -b $CHROMHMM -wa -wb |\
 awk -v COLS=$(head -1 $FILE | awk '{ print NF }') '{ tot[$(4+COLS)]+=1 } END { for (i in tot) print i"\t"tot[i] }' |\
 sort