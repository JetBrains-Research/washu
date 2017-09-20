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
    echo "Need 2 parameters! <BAM> <INSERT_SIZE>"
    exit 1
fi

BAM=$1
INSERT_SIZE=$2
SHIFT=$(($INSERT_SIZE / 2))

TMP_DIR=~/tmp
mkdir -p "${TMP_DIR}"

bedtools bamtobed -i ${BAM} |\
    awk -v OFS='\t' -v S=${SHIFT} \
    '{if ($6 != "-") {print($1, $2+S, $2+S+1)} else {if ($3-S>=1) {print($1, $3-S, $3-S+1)}}}' |\
    sort -k1,1 -k3,3n -k2,2n -T ${TMP_DIR}
