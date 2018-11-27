#!/bin/bash
# This script is used to compute union of peaks for given list of files.
#
# What happens:
# - Two peaks aver overlapping if they share at least one nucleotide
# - Prints only peaks that exist at least in one file.
#       4th column indicates tracks, parents of each peak in union.
#
# author Oleg Shpynov (oleg.shpynov@jetbrains.com)

which bedtools &>/dev/null || { echo "ERROR: bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
>&2 echo "union: $@"

if [[ $# -eq 0 ]]; then
  echo "ERROR: Empty arguments list"
  exit 1
fi

# Use temp file since folder can be read-only
TMP=$(mktemp)
if [ -f ${TMP} ]; then
    rm $TMP
fi

# Check configuration
[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util.sh
export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p "${TMPDIR}"

N=1
for FILE in $@
do
    NAME=${FILE##.*/}
    awk -v OFS='\t' -v N=${N}_${NAME} '{print $1,$2,$3,N}' ${FILE} >> ${TMP}
    N=$((N+1))
done

SORTED=$(mktemp)
sort -k1,1 -k2,2n -T ${TMPDIR} ${TMP} > ${SORTED}
bedtools merge -i ${SORTED} -c 4 -o distinct -delim "|"

# Cleanup
rm ${TMP} ${SORTED}
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir