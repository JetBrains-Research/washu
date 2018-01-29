#!/bin/bash
# This script is used to compute jaccard index of 2 BED files
# author Oleg Shpynov (oleg.shpynov@jetbrains.com)

which bedtools &>/dev/null || { echo "ERROR: bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }

>&2 echo "jaccard $@"

SORTED=NO
MERGED=NO
HELP=NO
###########################################
POSITIONAL=()

if [[ $# -eq 0 ]]; then
  HELP=YES
fi

while [[ $# -gt 0 ]]; do
key="$1"
case $key in
    --help|-h)
    HELP=YES
    shift # past argument with no value
    ;;

    -s) # Already sorted
    SORTED=YES
    shift # past argument with no value
    ;;

    -m) # Already merged
    MERGED=YES
    shift # past argument with no value
    ;;

    *)
    POSITIONAL+=("$key")
    shift
    ;;
esac
done
set -- "${POSITIONAL[@]}"
###########################################

if [ "${HELP}" == "YES" ]; then
  echo "Calculate Jaccard Index for A.bed B.bed"
  echo ""
  echo "Usage: jaccard.sh [OPTIONS] A.bed B.bed"
  echo ""
  echo "Options:"
  echo "  -s              Files already sorted, skip resort step"
  echo "  -m              Files already merged, skip merge step"
  echo "  -h|--help       Show help"
  echo ""
  echo "With no arguments - show help"
  exit 1
fi

if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <A.bed> <B.bed>"
    exit 1
fi

# Use temp file since folder can be read-only
TMP=$(mktemp)
if [ -f ${TMP} ]; then
    rm $TMP
fi

# Check configuration
[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util/util.sh
export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p "${TMPDIR}"

BED1=$1
BED2=$2

# 1. input files may contain intersecting intervals (e.g. introns vs transcripts)
# merging is obligatory there so as get correct result
#
# 2. whether files are sorted or not isn't critical here, but for sorted files
# we can use additional options which can speedup calculations

if [ "${MERGED}" == "NO" ]; then
    # merge command requires pre-sorted files
    if [ "${SORTED}" == "NO" ]; then
        BED1_SORTED=${TMPDIR}/1.sorted.bed
        BED2_SORTED=${TMPDIR}/2.sorted.bed
        sort -k1,1 -k2,2n -T ${TMPDIR} $BED1 > $BED1_SORTED
        sort -k1,1 -k2,2n -T ${TMPDIR} $BED2 > $BED2_SORTED
        BED1=$BED1_SORTED
        BED2=$BED2_SORTED
    fi

    BED1_MERGED=${TMPDIR}/1.merged.bed
    BED2_MERGED=${TMPDIR}/2.merged.bed
    bedtools merge -i $BED1 > $BED1_MERGED
    bedtools merge -i $BED2 > $BED2_MERGED
    BED1=$BED1_MERGED
    BED2=$BED2_MERGED
fi

if [ "${SORTED}" == "NO" ]; then
    SORTED_OPT=""
else
    SORTED_OPT=" -sorted"
fi

INTERSECT=$(bedtools intersect$SORTED_OPT -a $BED1 -b $BED2 | awk 'BEGIN{L=0}; {L+=$3-$2}; END{print(L)}')
UNION=$(bash ${WASHU_ROOT}/bed/union.sh $BED1 $BED2 | awk 'BEGIN{L=0}; {L+=$3-$2}; END{print(L)}')

# Empty union results in 0
if [[ $UNION -eq "0" ]]; then
    echo 0
else
    echo "$(bc -l <<< "$INTERSECT / $UNION")" | sed 's#^\.#0.#g'
fi

# TMP dir cleanup:
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir