#!/bin/bash
# This script is used to compute number of reads in peaks.
# author oleg.shpynov@jetbrains.com

which bedtools &>/dev/null || { echo "bedtools not found!"; exit 1; }
>&2 echo "rip.sh: $@"

if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <bam> <peaks>"
    exit 1
fi
 
READS_BAM=$1
PEAKS_FILE=$2

READS_PREFIX=${READS_BAM%%.bam}
PILEUP_BED=${READS_PREFIX}_pileup.bed
RIP_FILE=${PEAKS_FILE}_rip.csv

echo "BAM_FILE: ${READS_BAM}"
echo "PILEUP_BED: ${PILEUP_BED}"
echo "RIP_FILE: ${RIP_FILE}"

# If we already have rip file, do not recalculate
if [[ -f ${RIP_FILE} ]]; then
    echo "${RIP_FILE}"
    cat ${RIP_FILE}
    exit 0
fi

PEAKS_FILE_SORTED=${PEAKS_FILE}.sorted
INTERSECT_BED=${PEAKS_FILE}.intersect.bed

source $(dirname $0)/../parallel/util.sh 2> /dev/null
export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p "${TMPDIR}"

# To pileup bed
bedtools bamtobed -i ${READS_BAM} | awk -v OFS='\t' '{print $1,$2,$3}' > sort -k1,1 -k2,2n -T ${TMPDIR} -o ${PILEUP_BED}
READS=$(wc -l ${PILEUP_BED} | awk '{print $1}')

# To sorted bed
cat ${PEAKS_FILE} | awk -v OFS='\t' '{print $1,$2,$3}' | sort -k1,1 -k2,2n -T ${TMPDIR} > ${PEAKS_FILE_SORTED}
PEAKS=$(wc -l ${PEAKS_FILE_SORTED} | awk '{print $1}')

# Compute number of reads, intersecting with peaks
intersectBed -a ${PILEUP_BED} -b ${PEAKS_FILE_SORTED} -c -f 0.20 > ${INTERSECT_BED}
RIP=$(awk '{sum += $4} END {print sum}' ${INTERSECT_BED})

# Show result
echo "file,peaks_file,reads,peaks,rip" > ${RIP_FILE}
echo "${READS_BAM},${PEAKS_FILE},${READS},${PEAKS},${RIP}" >> ${RIP_FILE}

echo "${RIP_FILE}"
cat ${RIP_FILE}

# Cleanup
rm ${PILEUP_BED} ${PEAKS_FILE_SORTED} ${INTERSECT_BED}
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir