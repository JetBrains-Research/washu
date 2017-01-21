#!/bin/bash
# This script is used to compute overlap of peaks for given list of files.
#
# What happens:
# - Filter out unknown contigs and Y chromosome peaks
# - Two peaks aver overlapping if they share at least one nucleotide
# - For each of the combination of input files, number of overlap peaks are computed
#
# Example:
# > bash metapeaks.sh A.bed B.bed
# Output:
# Metapeaks A.bed B.bed
# PEAKS:       10       10
# 0 1	1
# 1 0	1
# 1 1	9
# author Konstantin Zaytsev
# author Oleg Shpynov

which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
>&2 echo "metapeaks: $@"

# FILTERED data on chromosomes only, i.e. no contig
CHRFILES=()
PEAKS=()
for i in $@
do
    tmpfile=${i}.chr_only.tmp
    grep -E "chr[0-9]+|chrX" $i | sort -k1,1 -k2,2n -k3,3n > $tmpfile
    CHRFILES+=("$tmpfile")
    peak=$(cat $tmpfile | wc -l)
    PEAKS+=("$peak")
done

echo "PEAKS: ${PEAKS[@]}"

range=$(seq -s, 6 1 $(($# + 5)))

METAPEAKS=$(
    multiIntersectBed -i "${CHRFILES[@]}" |\
    bedtools merge -c $range -o max |\
    # Zero problem: max of '0' is 2.225073859e-308 - known floating point issue in bedtools merge
    awk '{if (NR > 1) printf("\n"); printf("%s\t%s\t%s", $1, $2, $3); for (i=4; i<=NF; i++) printf("\t%d", int($i)); }' |\
    # Extract columns 4 up to the end
    awk '{for (i=4; i<=NF; i++) printf("%s%s", $i, (i==NF) ? "\n" : OFS)}' |\
    # Compute all the different lines and log it
    awk '{ tot[$0]++ } END { for (i in tot) print i"\t"tot[i] }' |\
    sort
)
echo "$METAPEAKS"
echo "$METAPEAKS" |\
# Short version for Venn Diagrams visualization - join last column with commas
awk 'NR > 1 { printf(", ") }{printf("%s", $NF)}'
echo

# Cleanup
rm ${CHRFILES[@]}