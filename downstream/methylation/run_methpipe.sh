#!/bin/bash

DESIGN=$1
TABLE=$2
OUTPUT=$3
NAME=$4

#rm ${NAME}.log

WDIR=`pwd`
ADJ_OUTPUT="${OUTPUT%.*}".narrow.adjusted.bed
REG_OUTPUT="${ADJ_OUTPUT%.*}".regions.bed

cd $WDIR
radmeth regression -factor case -o $OUTPUT $DESIGN $TABLE
radmeth adjust -bins 1:50:1 $OUTPUT > $ADJ_OUTPUT
radmeth merge -p 0.05 $ADJ_OUTPUT > $REG_OUTPUT