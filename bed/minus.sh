#!/bin/bash
# This script is used to compute *minus* of peaks for given 2 files.
#
# What happens:
# - Filter out unknown contigs
# - Two peaks aver overlapping if they share at least one nucleotide
# - Prints only peaks that is unique to first file
#
# author Oleg Shpynov (oleg.shpynov@jetbrains.com)

which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
>&2 echo "minus: $@"

if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <FILE1> <FILE2>"
    exit 1
fi

FILE1=$1
FILE2=$2

# FILTERED data on chromosomes only, i.e. no contig
grep -E "chr[0-9]+|chrX|chrY" $FILE1 | sort -k1,1 -k2,2n > ${FILE1}.tmp
grep -E "chr[0-9]+|chrX|chrY" $FILE2 | sort -k1,1 -k2,2n > ${FILE2}.tmp

multiIntersectBed -i ${FILE1}.tmp ${FILE2}.tmp |\
 bedtools merge -c 6,7 -o max |\
 # Zero problem: max of '0' is 2.225073859e-308 - known floating point issue in bedtools merge
 awk '{if (NR > 1) printf("\n"); printf("%s\t%s\t%s", $1, $2, $3); for (i=4; i<=NF; i++) printf("\t%d", int($i)); }' |\
 # NOTE[shpynov] use awk instead of grep, because grep has some problems with tab characters.
 awk "/\t1\t0/" |\
 awk '{for (i=1; i<=3; i++) printf("%s%s", $i, (i==3) ? "\n" : "\t")}' |\
 sort -k1,1 -k2,2n

# Cleanup
rm ${FILE1}.tmp ${FILE2}.tmp
