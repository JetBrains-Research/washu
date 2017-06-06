#!/usr/bin/env bash
# This script is used to select top N peaks for given MACS2 broad peaks
# Usage:
#   bash top_peaks.sh <PEAKS_FOLDER> <OUTPUT_FOLDER> <N>
#
# author Oleg Shpynov (oleg.shpynov@jetbrains.com)

>&2 echo "macs2_top_peaks.sh: $@"

if [ $# -lt 3 ]; then
    echo "Need 3 parameters! <PEAKS_FOLDER> <OUTPUT_FOLDER> <N>"
    exit 1
fi

PEAKS_FOLDER=$1
OUTPUT_FOLDER=$2
N=$3

if [[ ! -d ${PEAKS_FOLDER} ]]; then
    echo "Missing folder ${PEAKS_FOLDER}"
    exit 1
fi
mkdir -p ${OUTPUT_FOLDER}

for F in $(find ${PEAKS_FOLDER} -name '*.xls'| sed 's#./##g' | grep -v input); do
    echo $F
    cat $F | grep -E '#|name' > ${OUTPUT_FOLDER}/$F
    cat $F | grep -v -E '#|name' | sort -k8,8gr | head -$N >> ${OUTPUT_FOLDER}/$F
done
