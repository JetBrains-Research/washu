#!/usr/bin/env bash

if [ $# -lt 1 ]; then
    echo "Need 1 parameter! <folder_path>"
    echo ""
    echo "For given folder with peaks for each region 'r' from multi intersection calculates: abs(#{OD_i intersecting 'r'} - #{YD_i intersecting 'r'})."
    exit 1
fi

FOLDER=$1
cd ${FOLDER}

OD_FILES="$(find . -maxdepth 1 \( -wholename "*OD*island.bed" -or -wholename "*OD*Peak" -or \
    -wholename "*OD*_peaks.bed" \) | grep -v outlier)"

YD_FILES="$(find . -maxdepth 1 \( -wholename "*YD*island.bed" -or -wholename "*YD*Peak" -or \
    -wholename "*YD*_peaks.bed" \) | grep -v outlier)"

#echo "OD: "${OD_FILES}
#echo "YD: "${YD_FILES}
N_OD="$(echo "${OD_FILES}" | wc -l)"
echo "chr"$'\t'"start"$'\t'"end"$'\t'"cons"$'\t'"od_cons"$'\t'"yd_cons"$'\t'"absdiff"
bedtools multiinter -v -i ${OD_FILES} ${YD_FILES} | awk -v K="${N_OD}" 'function abs(v)
{return v < 0 ? -v : v}{
 ods=0; yds=0; for (i=6; i<6+K; i++) ods+=$i;
 for (i=6+K; i<=NF; i++) yds+=$i;
 print $1, $2, $3, $4, ods, yds, abs(ods-yds) }'
