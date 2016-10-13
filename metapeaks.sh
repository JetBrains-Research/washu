#!/bin/bash
# This script is used to compute overlap of peaks for given list of files.
#
# What happens:
# - Filter out unknown contigs and Y chromosome peaks
# - Two peaks aver overlapping if they share at least one nucleotide
# - For each of the combination of input files, number of overlap peaks are computed
# Example:
# > bash metapeaks.sh A.bed B.bed
# Output:
# PEAKS:       10       10
# 0 1	1
# 1 0	1
# 1 1	9
# author Konstantin Zaytsev
# author Oleg Shpynov

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

echo "PEAKS: ${PEAKS[@]}"

range=$(seq -s, 6 1 $(($# + 5)))

multiIntersectBed -i "${CHRFILES[@]}" |\
bedtools merge -c $range -o max |\
# Extract columns 4 up to the end
awk '{for (i=4; i<=NF; i++) printf("%s%s", $i, (i==NF) ? "\n" : OFS)}' |\
# Compute all the different lines and log it
awk '{ tot[$0]++ } END { for (i in tot) print i"\t"tot[i] }' |\
sort

# Cleanup
rm ${CHRFILES[@]}