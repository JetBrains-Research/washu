#!/bin/bash

# this script is from Tao Liu https://gist.github.com/taoliu/2469050 

which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
which bedGraphToBigWig &>/dev/null || { echo "bedGraphToBigWig not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }
which bedClip &>/dev/null || { echo "bedClip not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }
 
# end of checking
 
if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <bedgraph> <chrom info>"
    exit 1
fi
 
BDG_FILE=$1
CHROM_SIZES=$2

# Remove coordinates outside chromosome sizes
bedtools slop -i ${BDG_FILE} -g ${CHROM_SIZES} -b 0 | bedClip stdin ${CHROM_SIZES} ${BDG_FILE}.clip
# Fix problem with not sorted clip file
LC_COLLATE=C sort -k1,1 -k2,2n ${BDG_FILE}.clip > ${BDG_FILE}.sort.clip
bedGraphToBigWig ${BDG_FILE}.sort.clip ${CHROM_SIZES} ${BDG_FILE/bdg/bw}
# Cleanup
rm -f ${BDG_FILE}*.clip

