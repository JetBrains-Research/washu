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
    echo "Need 2 parameters! <BED_FILE> <CHROMHMM_MARKUP.bigBed | CHROMHMM_MARKUP.bed>"
    exit 1
fi
>&2 echo "chromhmm_stat $@"

FILE=$1
CHROMHMM=$2

# Convert to bed if required
if [[ ! $CHROMHMM == *.bed ]]; then
    NAME=${CHROMHMM%%.*}
    >&2 echo "ChromHMM markup is in bigBed format, converting to $CHROMHMM.bed"
    if [ ! -f "$CHROMHMM.bed" ]; then
        which bigBedToBed &>/dev/null || { echo "bigBedToBed not found! Download: http://hgdownload.cse.ucsc.edu/admin/exe/ or 'conda install ucsc-bigbedtobed'"; exit 1; }
        bigBedToBed $FILE $CHROMHMM.bed
    fi
    CHROMHMM=$CHROMHMM.bed
fi

# Compute intersection by chromhmm state
COLS=$(cat $FILE | grep "chr" | head -1 | awk '{ print NF }')
bedtools intersect -a "$FILE" -b "$CHROMHMM" -wa -wb |\
 awk -v COLS=$COLS -v OFS='\t' '{ states[$(4+COLS)]+=1 } END { for (i in states) print i,states[i] }' |\
 sort