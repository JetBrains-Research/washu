#!/bin/bash
# This script is used to compute number of reads in peaks.
# author oleg.shpynov@jetbrains.com

which bedtools &>/dev/null || { echo "bedtools not found!"; exit 1; }
>&2 echo "rip.sh: $@"

if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <bam> <peaks>"
    exit 1
fi
 
FILE=$1
NAME=${FILE%%.bam}
PEAKS=$2
 
# Compute number pileup bed
if [ ! -f "${NAME}_pileup.bed" ]; then
    # To pileup bed
    bedtools bamtobed -i ${FILE} > ${NAME}_pileup.bed
fi
READS=$(wc -l ${NAME}_pileup.bed | awk '{print $1}')

# Compute number of reads, intersecting with peaks
intersectBed -a ${NAME}_pileup.bed -b ${PEAKS} -c -f 0.20 > ${PEAKS}.intersectBed
RIP=$(awk '{sum += $7} END {print sum}' ${PEAKS}.intersectBed)

# Show result
echo "${FILE},${READS}" > ${PEAKS}_rip.txt
echo "${PEAKS},${RIP}" >> ${PEAKS}_rip.txt
cat ${PEAKS}_rip.txt

# Cleanup
rm ${NAME}_pileup.bed ${PEAKS}.intersectBed