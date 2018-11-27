#!/bin/bash
# This is a part of ENCODE-DCC chip-seq-pipeline, distributed with MIT license.
# https://github.com/ENCODE-DCC/chip-seq-pipeline/blob/master/dnanexus/filter_qc/src/filter_qc.py
# author oleg.shpynov@jetbrains.com

which bedtools &>/dev/null || { echo "ERROR: bedtools not found!"; exit 1; }
>&2 echo "pbc_nrf.sh: $@"

if [ $# -lt 1 ]; then
    echo "Need 1 parameter! <bam>"
    exit 1
fi
 
BAM=$1
>&2 echo "BAM: ${BAM}"

[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util.sh
export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p "${TMPDIR}"

PILEUP_BED=$(pileup ${BAM})
>&2 echo "PILEUP_BED: ${PILEUP_BED}"

T=$'\t'
>&2 echo "TotalReadPairs${T}DistinctReadPairs${T}OneReadPair${T}TwoReadPairs${T}\
NRF=Distinct/Total${T}PBC1=OnePair/Distinct${T}PBC2=OnePair/TwoPair"

cat ${PILEUP_BED} | \
    sort -k1,1 -k3,3n -k2,2n -k6,6 -T ${TMPDIR} | \
    awk -v OFS='\t' '{print $1,$2,$3,$6}' | uniq -c | \
    awk 'BEGIN{mt=0;m0=0;m1=0;m2=0}
    ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1}
    END{
        if (mt!=0){m0_t=m0/mt} else {m0_t=-1.0};
        if (m0!=0){m1_0=m1/m0} else {m1_0=-1.0};
        if (m2!=0){m1_2=m1/m2} else {m1_2=-1.0};
        printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0_t,m1_0,m1_2;
    }'

# Cleanup
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir