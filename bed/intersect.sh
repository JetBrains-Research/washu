#!/bin/bash
# This script is used to compute intersection of peaks for given list of files.
#
# What happens:
# - Filter out unknown contigs and Y chromosome peaks
# - Two peaks aver overlapping if they share at least one nucleotide
# - Prints only peaks that overlap in all files (merged)
#
# author Oleg Shpynov (oleg.shpynov@jetbrains.com)

which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
>&2 echo "intersect: $@"

# FILTERED data on chromosomes only, i.e. no contig
CHRFILES=()
PEAKS=()
for i in $@;
do
    tmpfile=${i}.chr_only.tmp
    grep -E "chr[0-9]+|chrX" $i > $tmpfile
    CHRFILES+=("$tmpfile")
    peak=$(cat $tmpfile | wc -l)
    PEAKS+=("$peak")
done

range=$(seq -s, 6 1 $(($# + 5)))
pattern=$(printf '\t1%.0s' $(seq 1 $#))

multiIntersectBed -i "${CHRFILES[@]}" |\
 bedtools merge -c $range -o max |\
 # Zero problem: max of '0' is 2.225073859e-308 - known floating point issue in bedtools merge
 awk '{if (NR > 1) printf("\n"); printf("%s\t%s\t%s", $1, $2, $3); for (i=4; i<=NF; i++) printf("\t%d", int($i)); }' |\
 # NOTE[shpynov] use awk instead of grep, because grep has some problems with tab characters.
 awk "/$pattern/" |\
 awk '{for (i=1; i<=3; i++) printf("%s%s", $i, (i==3) ? "\n" : "\t")}' |\
 sort -k1,1 -k2,2n

# Cleanup
rm ${CHRFILES[@]}
