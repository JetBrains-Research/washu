#!/bin/bash
# This is a part of ENCODE-DCC chip-seq-pipeline, distributed with MIT license.
# https://github.com/ENCODE-DCC/chip-seq-pipeline/blob/master/dnanexus/filter_qc/src/filter_qc.py
# author oleg.shpynov@jetbrains.com

which bedtools &>/dev/null || { echo "bedtools not found!"; exit 1; }
>&2 echo "pbc_nrf.sh: $@"

if [ $# -lt 1 ]; then
    echo "Need 1 parameter! <bam>"
    exit 1
fi
 
BAM=$1
PILEUP_BED=${BAM/.bam/_pileup.bed}

>&2 echo "BAM: ${BAM}"
>&2 echo "PILEUP_BED: ${PILEUP_BED}"

source $(dirname $0)/../parallel/util.sh 2> /dev/null
export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p "${TMPDIR}"

# To pileup bed
if [ -f ${PILEUP_BED} ]; then
    >&2 echo "Pileup file already exists: ${PILEUP_BED}"
else
    # Safely recalc in tmpfile if not exists:
    # - recalc in tmp file
    # - then move to desired location
    # This script could be launched in parallel on cluster, and *pileup.bed fill
    # is shared among peaks filtering tasks. So as not to get inconsistent
    # file state let's change file in atomic like way
    PILEUP_TMP=$(mktemp $TMPDIR/pileup.XXXXXX.bed)
    >&2 echo "Calculate pileup file in tmp file: ${PILEUP_TMP}"
    bedtools bamtobed -i ${BAM} > ${PILEUP_TMP}
    if [ -f ${PILEUP_BED} ]; then
        >&2 echo "  Ignore result, file has been already calculated: ${PILEUP_BED}"
    else
        # if still doesn't exists:
        >&2 echo "  Move ${PILEUP_TMP} -> ${PILEUP_BED}"
        mv ${PILEUP_TMP} ${PILEUP_BED}
    fi
fi

# PBC File output:
#   TotalReadPairs <tab>
#   DistinctReadPairs <tab>
#   OneReadPair <tab>
#   TwoReadPairs <tab>
#   NRF=Distinct/Total <tab>
#   PBC1=OnePair/Distinct <tab>
#   PBC2=OnePair/TwoPair
cat ${PILEUP_BED} | grep -E 'chr[0-9XY]+' | \
    awk  -v OFS='\t' '{print $1,$2,$3,$6}' | \
    sort | uniq -c | \
    awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}'

# Cleanup
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir