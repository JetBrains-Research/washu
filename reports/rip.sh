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
PILEUP_BED=${NAME}_pileup.bed

PEAKS_FILE=$2
RIP_FILE=${PEAKS_FILE}_rip.csv
PEAKS_FILE_SORTED=${PEAKS_FILE}.sorted
INTERSECT_BED=${PEAKS_FILE}.intersect.bed

# To pileup bed
bedtools bamtobed -i ${FILE} | awk -v OFS='\t' '{print $1,$2,$3}' | sort -k1,1 -k2,2n > ${PILEUP_BED}
READS=$(wc -l ${PILEUP_BED} | awk '{print $1}')

# To sorted bed
cat ${PEAKS_FILE} | awk -v OFS='\t' '{print $1,$2,$3}' | sort -k1,1 -k2,2n > ${PEAKS_FILE_SORTED}
PEAKS=$(wc -l ${PEAKS_FILE_SORTED} | awk '{print $1}')

# Compute number of reads, intersecting with peaks
intersectBed -a ${PILEUP_BED} -b ${PEAKS_FILE_SORTED} -c -f 0.20 > ${INTERSECT_BED}
RIP=$(awk '{sum += $4} END {print sum}' ${INTERSECT_BED})

# Show result
echo "file,peaks_file,reads,peaks,rip" > ${RIP_FILE}
echo "${FILE},${PEAKS_FILE},${READS},${PEAKS},${RIP}" >> ${RIP_FILE}

echo "${RIP_FILE}"
cat ${RIP_FILE}

# Cleanup
rm ${PILEUP_BED} ${PEAKS_FILE_SORTED} ${INTERSECT_BED}