#!/usr/bin/env bash
# author: oleg.shpynov@jetbrains.com

# Check tool.
which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }

if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <BED_FILE> <GENES.ANNOTATION.gtf | GENES.ANNOTATION.bed>"
    echo "Download annotation at: https://www.gencodegenes.org/"
    exit 1
fi
>&2 echo "closest_gene $@"

FILE=$1
GENES=$2

if [[ ! $GENES == *.bed ]]; then
    # Gtf to sorted tsv conversion
    GENES_TSV=${GENES/gtf/sorted.tsv}
    if [ ! -f ${GENES_TSV} ]; then
        >&2 echo "Converting gtf to ${GENES_TSV}"
        GENE_NAME_FIELD=$(cat ${GENES} | grep "chr1" | head -1 | awk '{for (i=1; i<NF; i++) {if ($i=="gene_name") print (i+1)}}')
        cat ${GENES} |  awk -v GN=${GENE_NAME_FIELD} 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$GN}}' | tr -d '";' |\
         sort -k1,1 -k2,2n > ${GENES_TSV}
    fi
    GENES=${GENES_TSV}
fi

COLS=$(cat $FILE | grep "chr" | head -1 | awk '{ print NF }')
bedtools closest -a ${FILE} -b ${GENES} -d | awk -v COLS=$COLS -v OFS='\t' '{ print $1,$2,$3,$(COLS+4),$(COLS+5) }' | sort