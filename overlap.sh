#!/bin/bash
# This script is used to compute overlap of peaks for given list of files.
#
# What happens:
# - Filter out unknown contigs and Y chromosome peaks
# - Two peaks aver overlapping if they share at least one nucleotide
# - Prints only peaks that overlap in all files (merged)
#
# author Oleg Shpynov (oleg.shpynov@jetbrains.com)

which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }

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

keys=$(for i in $(seq 4 $(($# + 3)));
do
    echo "-k$i,$i"
done
)
range=$(seq -s, 6 1 $(($# + 5)))

multiIntersectBed -i "${CHRFILES[@]}" |\
bedtools merge -c $range -o max |\
awk -v OFS="\t" '{for (i=4; i<= NF; i++) $i = int($i); print $0}' |\
sort $keys |\
grep -e "\t1\t1" |\
awk -v OFS="\t" '{for (i=1; i<=3; i++) printf("%s%s", $i, (i==3) ? "\n" : OFS)}'

# Cleanup
rm $CHRFILES
