#!/usr/bin/env bash
# author: oleg.shpynov@jetbrains.com

# Check tools
which bedtools &>/dev/null || {
    echo "bedtools not found! You can install it using:"
    echo "  conda install -c bioconda bedtools"
    echo "For further details see http://code.google.com/p/bedtools"
    exit 1;
   }

>&2 echo "bam2reads $@"
if [ $# -lt 1 ]; then
    echo "Need 1 parameters! <BAM>"
    exit 1
fi

BAM=$1

bedtools bamtobed -i $BAM | awk -v OFS='\t' '{print $1,$2,$3,$6}'
