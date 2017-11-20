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

PEAKS_FILE_SORTED=${PEAKS_FILE}.sorted
INTERSECT_BED=${PEAKS_FILE}.intersect.bed

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
    TMPFILE=$(mktemp $TMPDIR/pileup.XXXXXX.bed)
    echo "Calculate pileup file in tmp file: ${TMPFILE}"
    bedtools bamtobed -i ${READS_BAM} | awk -v OFS='\t' '{print $1,$2,$3}' | sort -k1,1 -k2,2n -T ${TMPDIR} -o ${TMPFILE}
    if [ -f ${PILEUP_BED} ]; then
        echo "  Ignore result, file has been already calculated by smb else: ${PILEUP_BED}"
    else
        # if still doesn't exists:
        echo "  Move ${TMPFILE} -> ${PILEUP_BED}"
        mv ${TMPFILE} ${PILEUP_BED}
    fi
fi
READS=$(wc -l ${PILEUP_BED} | awk '{print $1}')
echo "READS: $READS"

# To sorted bed
cat ${PEAKS_FILE} | awk -v OFS='\t' '{print $1,$2,$3}' | sort -k1,1 -k2,2n -T ${TMPDIR} > ${PEAKS_FILE_SORTED}
PEAKS=$(wc -l ${PEAKS_FILE_SORTED} | awk '{print $1}')
echo "PEAKS: $PEAKS"

# Compute number of reads, intersecting with peaks
intersectBed -a ${PILEUP_BED} -b ${PEAKS_FILE_SORTED} -c -f 0.20 > ${INTERSECT_BED}
# _pileup.bed can have different number of columns
COLS=$(cat ${PILEUP_BED} | head -1 | awk '{ print NF }')
RIP=$(awk -v COLS=$COLS '{sum += $(COLS+1)} END {print sum}' ${INTERSECT_BED})
echo "RIP: $RIP"

echo "file,peaks_file,reads,peaks,rip" > ${RIP_FILE}
echo "${READS_BAM},${PEAKS_FILE},${READS},${PEAKS},${RIP}" >> ${RIP_FILE}

# Cleanup
rm ${PEAKS_FILE_SORTED} ${INTERSECT_BED}
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir