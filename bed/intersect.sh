#!/bin/bash
# This script is used to compute intersection of peaks for given list of files.
#
# What happens:
# - Two peaks aver overlapping if they share at least one nucleotide
# - Prints only peaks that overlap in all files (merged)
#
# author Oleg Shpynov (oleg.shpynov@jetbrains.com)

which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
>&2 echo "intersect: $@"

if [[ $# -eq 0 ]]; then
  echo "ERROR: Empty arguments list"
  exit 1
fi

# Optional load technical stuff:
source $(dirname $0)/../parallel/util/util.sh 2> /dev/null
export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p "${TMPDIR}"

SORTED_FILES=()
for F in $@
do
    # Folder with source file be read-only, use temp file
    SORTED=$(mktemp)
    sort -k1,1 -k2,2n -T ${TMPDIR} $F > ${SORTED}
    SORTED_FILES+=("$SORTED")
done

range=$(seq -s, 6 1 $(($# + 5)))
pattern=$(printf '\t1%.0s' $(seq 1 $#))

bedtools multiinter -i "${SORTED_FILES[@]}" |\
 bedtools merge -c $range -o max |\
 # Zero problem: max of '0' is 2.225073859e-308 - known floating point issue in bedtools merge
 awk '{if (NR > 1) printf("\n"); printf("%s\t%s\t%s", $1, $2, $3); for (i=4; i<=NF; i++) printf("\t%d", int($i)); }' |\
 # NOTE[shpynov] use awk instead of grep, because grep has some problems with tab characters.
 awk "/$pattern$/" |\
 awk '{for (i=1; i<=3; i++) printf("%s%s", $i, (i==3) ? "\n" : "\t")}' |\
 sort -k1,1 -k2,2n -T ${TMPDIR}

# Cleanup
rm ${SORTED_FILES[@]}
# TMP dir cleanup:
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir

