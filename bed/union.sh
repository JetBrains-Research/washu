#!/bin/bash
# This script is used to compute union of peaks for given list of files.
#
# What happens:
# - Two peaks aver overlapping if they share at least one nucleotide
# - Prints only peaks that exist at least in one file.
#       4th column indicates tracks, parents of each peak in union.
#
# author Oleg Shpynov (oleg.shpynov@jetbrains.com)

which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
>&2 echo "union: $@"

# Use temp file since folder can be read-only
TMP=$(mktemp)
if [ -f ${TMP} ]; then
    rm $TMP
fi

N=1
for FILE in $@
do
    NAME=${FILE##.*/}
    awk -v OFS='\t' -v N=${N}_${NAME} '{print $1,$2,$3,N}' $FILE >> ${TMP}
    N=$((N+1))
done

SORTED=$(mktemp)
sort -k1,1 -k2,2n ${TMP} > ${SORTED}
bedtools merge -i ${SORTED} -c 4 -o collapse -delim "|"

# Cleanup
rm ${TMP} ${SORTED}
