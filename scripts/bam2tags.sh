#!/usr/bin/env bash
# author: oleg.shpynov@jetbrains.com

# Check tools
which bedtools &>/dev/null || {
    echo "bedtools not found! You can install it using:"
    echo "  conda install -c bioconda bedtools"
    echo "For further details see http://code.google.com/p/bedtools"
    exit 1;
   }

>&2 echo "bam2tags $@"
if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <BAM> <FRAGMENT>"
    exit 1
fi

BAM=$1
FRAGMENT=$2
SHIFT=$(($FRAGMENT / 2))

# Check configuration
[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util/util.sh
export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p "${TMPDIR}"

PILEUP_BED=$(pileup ${BAM})
cat ${PILEUP_BED} |\
    awk -v OFS='\t' -v S=${SHIFT} \
    '{if ($6 != "-") {print($1, $2+S, $2+S+1)} else {if ($3-S>=1) {print($1, $3-S, $3-S+1)}}}' |\
    sort -u -k1,1 -k3,3n -k2,2n -T ${TMPDIR}

# TMP dir cleanup:
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir