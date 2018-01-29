#!/bin/bash
# This script is used to compute overlap of peaks for given list of files.
#
# What happens:
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

which bedtools &>/dev/null || { echo "ERROR: bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
>&2 echo "metapeaks: $@"

if [[ $# -eq 0 ]]; then
  echo "ERROR: Empty arguments list"
  exit 1
fi

# Check configuration
[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util/util.sh
export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p "${TMPDIR}"

SORTED_FILES=()
PEAKS=""
T=$'\t'
for F in $@
do
    # sometimes we need to trim 'wc -l' results, at least on my mac
    # the easiest suitable way here is using 'xargs'
    # Read about prons and cons at https://stackoverflow.com/questions/369758/how-to-trim-whitespace-from-a-bash-variable
    NPEAKS="$(cat ${F} | wc -l | xargs echo)"
    PEAKS="$PEAKS$T$NPEAKS"
    # Folder with source file be read-only, use temp file
    SORTED=$(mktemp)
    sort -k1,1 -k2,2n -T ${TMPDIR} $F > ${SORTED}
    SORTED_FILES+=("$SORTED")
done

echo "PEAKS:${PEAKS[@]}"

range=$(seq -s, 6 1 $(($# + 5)))

METAPEAKS=$(
    bedtools multiinter -i "${SORTED_FILES[@]}" |\
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
rm ${SORTED_FILES[@]}
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir
