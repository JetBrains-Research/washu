#!/bin/bash
# This script is used to compute jaccard index of 2 BED files
# author Oleg Shpynov (oleg.shpynov@jetbrains.com)

which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }

>&2 echo "jaccard $@"
if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <BED1> <BED2>"
    exit 1
fi

# Use temp file since folder can be read-only
TMP=$(mktemp)
if [ -f ${TMP} ]; then
    rm $TMP
fi

# Optional load technical stuff:
source $(dirname $0)/../parallel/util.sh 2> /dev/null
export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p "${TMPDIR}"

BED1=${TMPDIR}/1.bed
sort -k1,1 -k2,2n -T ${TMPDIR} $1 > $BED1
BED2=${TMPDIR}/2.bed
sort -k1,1 -k2,2n -T ${TMPDIR} $2 > $BED2

INTERSECT=$(bedtools intersect -a $BED1 -b $BED2 | awk 'BEGIN{L=0}; {L+=$3-$2}; END{print(L)}')
UNION=$(bash $(dirname $0)/union.sh $BED1 $BED2 | awk 'BEGIN{L=0}; {L+=$3-$2}; END{print(L)}')
echo "0$(bc -l <<< "$INTERSECT / $UNION")"


# TMP dir cleanup:
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir