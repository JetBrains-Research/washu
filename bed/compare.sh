#!/bin/bash
# This script is used to compute intersection of peaks for given list of files.
#
# What happens:
# - Filter out unknown contigs and Y chromosome peaks
# - Two peaks aver overlapping if they share at least one nucleotide
# - Produces 3 files: exclusive peaks for condition 1 and 2 and common peaks.
#
# author Oleg Shpynov (oleg.shpynov@jetbrains.com)

which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }

if [ $# -lt 3 ]; then
    echo "Need 3 parameters! <BED_1> <BED_2> <OUT_PREFIX>"
    exit 1
fi
>&2 echo "compare $@"

BED_1=$1
BED_2=$2
OUT_PREFIX=$3

# FILTERED data on chromosomes only, i.e. no contig
FILES=( $BED_1 $BED_2 )
CHRFILES=()
for i in ${FILES[@]}
do
    tmpfile=${i}.chr_only.tmp
    grep -E "chr[0-9]+|chrX" $i | sort -k1,1 -k2,2n > $tmpfile
    CHRFILES+=("$tmpfile")
done

TMP_FILE=${OUT_PREFIX}_all.txt

# Compute common and exclusive peaks
multiIntersectBed -i ${CHRFILES[@]} |\
bedtools merge -c 6,7 -o max |\
# Zero problem: max of '0' is 2.225073859e-308 - known floating point issue in bedtools merge
 awk -v OFS="\t" '{print $1,$2,$3,int($4),int($5)}' > ${TMP_FILE}

cat ${TMP_FILE} | awk '/\t1\t0/' |\
 awk '{for (i=1; i<=3; i++) printf("%s%s", $i, (i==3) ? "\n" : "\t")}' |\
 sort -k1,1 -k2,2n > ${OUT_PREFIX}_cond1.bed

cat ${TMP_FILE} | awk '/\t0\t1/' |\
 awk '{for (i=1; i<=3; i++) printf("%s%s", $i, (i==3) ? "\n" : "\t")}' |\
 sort -k1,1 -k2,2n > ${OUT_PREFIX}_cond2.bed

cat ${TMP_FILE} | awk '/\t1\t1/' |\
 awk '{for (i=1; i<=3; i++) printf("%s%s", $i, (i==3) ? "\n" : "\t")}' |\
 sort -k1,1 -k2,2n > ${OUT_PREFIX}_common.bed

# Cleanup
rm ${TMP_FILE}
rm ${CHRFILES[@]}