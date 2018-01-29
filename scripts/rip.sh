#!/bin/bash
# This script is used to compute number of reads in peaks.
# author oleg.shpynov@jetbrains.com

which bedtools &>/dev/null || { echo "ERROR: bedtools not found!"; exit 1; }
>&2 echo "rip.sh: $@"

if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <bam> <peaks>"
    exit 1
fi
 
READS_BAM=$1
PEAKS_FILE=$2

READS_PREFIX=${READS_BAM%%.bam}
RIP_FILE=${PEAKS_FILE}_rip.csv

echo "BAM_FILE: ${READS_BAM}"
echo "PEAKS_FILE: ${PEAKS_FILE}"
echo "RIP_FILE: ${RIP_FILE}"

# If we already have rip file, do not recalculate
if [[ -f ${RIP_FILE} ]]; then
    echo "${RIP_FILE}"
    cat ${RIP_FILE}
    exit 0
fi

[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util/util.sh
export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p "${TMPDIR}"

PILEUP_BED=$(pileup ${READS_BAM})
echo "PILEUP_BED: ${PILEUP_BED}"

READS=$(wc -l ${PILEUP_BED} | awk '{print $1}')
echo "READS: $READS"

PEAKS=$(wc -l ${PEAKS_FILE} | awk '{print $1}')
echo "PEAKS: $PEAKS"

LENGTH=$(awk '{l+=($3-$2)} END {print l}' ${PEAKS_FILE})
echo "LENGTH: $LENGTH"

# Compute number of reads, intersecting with peaks
INTERSECT_TMP=$(mktemp $TMPDIR/intersect.XXXXXX.bed)
echo "Calculate intersection in tmp file: ${INTERSECT_TMP}"
bedtools intersect -a ${PILEUP_BED} -b ${PEAKS_FILE} -c -f 0.20 > ${INTERSECT_TMP}

# _pileup.bed can have different number of columns
COLS=$(cat ${PILEUP_BED} | head -n 1 | awk '{ print NF }')
RIP=$(awk -v COLS=$COLS '{sum += $(COLS+1)} END {print sum}' ${INTERSECT_TMP})
echo "RIP: $RIP"

echo "file,peaks_file,reads,peaks,length,rip" > ${RIP_FILE}
echo "${READS_BAM},${PEAKS_FILE},${READS},${PEAKS},${LENGTH},${RIP}" >> ${RIP_FILE}

# Cleanup
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir