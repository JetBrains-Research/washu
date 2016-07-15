#!/usr/bin/env bash

# Original code: http://crazyhottommy.blogspot.ru/2014/10/convert-bam-file-to-bigwig-file-and.html
# author: oleg.shpynov@jetbrains.com

# Check tool.
which genomeCoverageBed &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }

if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <BAM> <chrom.sizes>"
    exit 1
fi

BAM=$1
CHROM_SIZES=$2

if [ ! -f "${CHROM_SIZES}" ]; then
  echo "File not found: ${CHROM_SIZES}"
  exit 1
fi

NAME=${BAM%%.bam}
genomeCoverageBed -ibam $BAM -bg -g ${CHROM_SIZES} > ${NAME}.bdg
bash ~/work/washu/bdg2bw.sh ${NAME}.bdg ${CHROM_SIZES}

