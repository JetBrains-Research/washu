#!/bin/bash

# this script is from Tao Liu https://gist.github.com/taoliu/2469050 
# check commands: slopBed, bedGraphToBigWig and bedClip
 
which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
which bedGraphToBigWig &>/dev/null || { echo "bedGraphToBigWig not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }
which bedClip &>/dev/null || { echo "bedClip not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }
 
# end of checking
 
if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <bedgraph> <chrom info>"
    exit 1
fi
 
F=$1
G=$2

# Remove coordinates outside chromosome sizes
bedtools slop -i ${F} -g ${G} -b 0 | bedClip stdin ${G} ${F}.clip
# Fix problem with not sorted clip file
LC_COLLATE=C sort -k1,1 -k2,2n ${F}.clip > ${F}.sort.clip
bedGraphToBigWig ${F}.sort.clip ${G} ${F/bdg/bw}
 
rm -f ${F}*.clip

