#!/bin/bash
# This script is used to compute union of peaks for given list of files.
#
# What happens:
# - Filter out unknown contigs and Y chromosome peaks
# - Two peaks aver overlapping if they share at least one nucleotide
# - Prints only peaks that exist at least in one file.
#       4th column indicates tracks, parents of each peak in union.
#
# author Oleg Shpynov (oleg.shpynov@jetbrains.com)

which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
>&2 echo "union: $@"

TMP="_to_merge.bed"
if [ -f ${TMP} ]; then
    rm $TMP
fi

# FILTERED data on chromosomes only, i.e. no contig
n=1
for f in $@
do
    grep -E "chr[0-9]+|chrX" $f | awk -v OFS='\t' -v N=$n '{print $1,$2,$3,N}' >> ${TMP}
    n=$((n+1))
done
sort -k1,1 -k2,2n -k3,3n ${TMP} > ${TMP}.sorted
bedtools merge -i ${TMP}.sorted -c 4 -o collapse -delim "|"

# Cleanup
rm ${TMP} ${TMP}.sorted
