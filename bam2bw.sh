#!/usr/bin/env bash

# Original code: http://crazyhottommy.blogspot.ru/2014/10/convert-bam-file-to-bigwig-file-and.html
# author: oleg.shpynov@jetbrains.com

# Check tools
which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
which bedGraphToBigWig &>/dev/null || { echo "bedGraphToBigWig not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }
which bedClip &>/dev/null || { echo "bedClip not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }

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

bedtools genomecov -ibam $BAM -bg -g ${CHROM_SIZES} | bedClip stdin ${CHROM_SIZES} ${NAME}.bdg
# CLIP output should be sorted for bedGraphToBigWig call
LC_COLLATE=C sort -k1,1 -k2,2n ${NAME}.bdg > ${NAME}.sorted.bdg
bedGraphToBigWig ${NAME}.sorted.bdg ${CHROM_SIZES} ${NAME}.bw

# Cleanup
rm -f ${NAME}*.bdg