#!/usr/bin/env bash
# author: oleg.shpynov@jetbrains.com

# Check tool.
which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }

if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <BED_FILE> <GENES.ANNOTATION.gtf | GENES.ANNOTATION.bed>"
    exit 1
fi
>&2 echo "closest_gene $@"

FILE=$1
GENES=$2

if [[ ! $GENES == *.bed ]]; then
    # Gtf to sorted bed conversion
    GENES_BED=${GENES/gtf/sorted.bed}
    if [ ! -f ${GENES_BED} ]; then
        >&2 echo "Converting gtf to ${GENES_BED}"
        cat ${GENES} |  awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$10,".",$7}}' | tr -d '";' |\
         sort -k1,1 -k2,2n > ${GENES_BED}
    fi
    GENES=${GENES_BED}
fi

COLS=$(cat $FILE | grep "chr" | head -1 | awk '{ print NF }')
bedtools closest -a ${FILE} -b ${GENES} | awk -v COLS=$COLS '{ print $(4+COLS) }' | sort | uniq