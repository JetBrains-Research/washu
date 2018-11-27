#!/bin/bash
# This script is used to compute number of reads in peaks.
# author oleg.shpynov@jetbrains.com

which bedtools &>/dev/null || { echo "ERROR: bedtools not found!"; exit 1; }
>&2 echo "rip.sh: $@"

if [[ $# -lt 2 ]]; then
    echo "Need 2 parameters! <reads> <peaks>"
    exit 1
fi

READS=$1
PEAKS_FILE=$2

RIP_FILE=${PEAKS_FILE}_rip.csv

echo "READS_FILE: ${READS}"
echo "PEAKS_FILE: ${PEAKS_FILE}"
echo "RIP_FILE: ${RIP_FILE}"

# If we already have rip file, do not recalculate
if [[ -f ${RIP_FILE} ]]; then
    echo "${RIP_FILE}"
    cat ${RIP_FILE}
    exit 0
fi

[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util.sh
export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p "${TMPDIR}"

if (echo ${READS} | grep -q ".*\\.bam$"); then
    READS_PREFIX=${READS%%.bam}
    PILEUP_BED=$(pileup ${READS})
elif (echo ${READS} | grep -q ".*\\.bed$"); then
    READS_PREFIX=${READS%%.bed}
    PILEUP_BED=${READS}
elif (echo ${READS} | grep -q ".*\\.bed.gz$"); then
    READS_PREFIX=${READS%%.bed.gz}
    PILEUP_BED=${READS%%.gz}
    if [[ ! -e ${PILEUP_BED} ]]; then
        gunzip -c ${READS} > ${PILEUP_BED}
    fi;
else
    echo "Unrecognized format of $READS, expected .bam, .bed or .bed.gz"; exit 1;
fi;

echo "PILEUP_BED: ${PILEUP_BED}"

READS=$(wc -l ${PILEUP_BED} | awk '{print $1}')
echo "READS: $READS"

PEAKS=$(wc -l ${PEAKS_FILE} | awk '{print $1}')
echo "PEAKS: $PEAKS"

LENGTH=$(awk 'BEGIN{l=0} {l+=($3-$2)} END{print l}' ${PEAKS_FILE})
echo "LENGTH: $LENGTH"

# Compute number of reads, intersecting with peaks
INTERSECT_TMP=$(mktemp $TMPDIR/intersect.XXXXXX.bed)
echo "Calculate intersection in tmp file: ${INTERSECT_TMP}"
bedtools intersect -a ${PILEUP_BED} -b ${PEAKS_FILE} -c -f 0.20 > ${INTERSECT_TMP}

# _pileup.bed can have different number of columns
COLS=$(cat ${PILEUP_BED} | head -n 1 | awk '{ print NF }')
RIP=$(awk -v COLS=${COLS} '{sum += $(COLS+1)} END {print sum}' ${INTERSECT_TMP})
echo "RIP: $RIP"

echo "file,peaks_file,reads,peaks,length,rip" > ${RIP_FILE}
echo "${READS},${PEAKS_FILE},${READS},${PEAKS},${LENGTH},${RIP}" >> ${RIP_FILE}

# Cleanup
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir