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
echo "PEAKS_FILE: ${PEAKS_FILE}"
echo "PILEUP_BED: ${PILEUP_BED}"
echo "RIP_FILE: ${RIP_FILE}"

# If we already have rip file, do not recalculate
if [[ -f ${RIP_FILE} ]]; then
    echo "${RIP_FILE}"
    cat ${RIP_FILE}
    exit 0
fi

source $(dirname $0)/../parallel/util.sh 2> /dev/null
export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p "${TMPDIR}"

# To pileup bed
if [ -f ${PILEUP_BED} ]; then
    echo "Pileup file already exists: ${PILEUP_BED}"
else
    # Safely recalc in tmpfile if not exists:
    # - recalc in tmp file
    # - then move to desired location
    # This script could be launched in parallel on cluster, and *pileup.bed fill
    # is shared among peaks filtering tasks. So as not to get inconsistent
    # file state let's change file in atomic like way
    PILEUP_TMP=$(mktemp $TMPDIR/pileup.XXXXXX.bed)
    echo "Calculate pileup file in tmp file: ${PILEUP_TMP}"
    bedtools bamtobed -i ${READS_BAM} > ${PILEUP_TMP}
    if [ -f ${PILEUP_BED} ]; then
        echo "  Ignore result, file has been already calculated: ${PILEUP_BED}"
    else
        # if still doesn't exists:
        echo "  Move ${PILEUP_TMP} -> ${PILEUP_BED}"
        mv ${PILEUP_TMP} ${PILEUP_BED}
    fi
fi
READS=$(wc -l ${PILEUP_BED} | awk '{print $1}')
echo "READS: $READS"

PEAKS=$(wc -l ${PEAKS_FILE} | awk '{print $1}')
echo "PEAKS: $PEAKS"

# Compute number of reads, intersecting with peaks
INTERSECT_TMP=$(mktemp $TMPDIR/intersect.XXXXXX.bed)
echo "Calculate intersection in tmp file: ${INTERSECT_TMP}"
bedtools intersect -a ${PILEUP_BED} -b ${PEAKS_FILE} -c -f 0.20 > ${INTERSECT_TMP}

# _pileup.bed can have different number of columns
COLS=$(cat ${PILEUP_BED} | head -1 | awk '{ print NF }')
RIP=$(awk -v COLS=$COLS '{sum += $(COLS+1)} END {print sum}' ${INTERSECT_TMP})
echo "RIP: $RIP"

echo "file,peaks_file,reads,peaks,rip" > ${RIP_FILE}
echo "${READS_BAM},${PEAKS_FILE},${READS},${PEAKS},${RIP}" >> ${RIP_FILE}

# Cleanup
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir