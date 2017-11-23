#!/bin/bash
# This script is used to compute BED files consensus
# author Petr Tsurinov (petr.tsurinov@jetbrains.com)

which bedtools &>/dev/null || {
    echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>";
    exit 1;
    }

>&2 echo "consensus $@"

PERCENT=0
COUNT=0
HELP=NO
TOOL="undefined"
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

    -p) # Consensus will be taken as percent
    PERCENT="$2"
    shift 2 # past argument
    ;;

    -c) # Consensus will be taken as tracks count
    COUNT="$2"
    shift 2 # past argument
    ;;

    *)
    POSITIONAL+=("$key")
    shift
    ;;
esac
done
set -- "${POSITIONAL[@]}"
###########################################

if [ "${HELP}" == "YES" ] || [ ${PERCENT} -lt 0 ] || [ ${COUNT} -lt 0 ]
    ([ ${PERCENT} -le 0 ] && [ ${COUNT} -le 0 ]) ||
    ([ ${PERCENT} -gt 0 ] && [ ${COUNT} -gt 0 ]); then
  echo "Calculate consensus for peaks in selected folder"
  echo ""
  echo "Usage: consensus.sh [OPTIONS] folder_path"
  echo ""
  echo "Options:"
  echo "  -p number       Consensus percent should be taken as number (cannot be used with -c)"
  echo "  -c number       Count of tracks for consensus should be taken as number
                          (cannot be used with -p)"
  echo "  -h|--help       Show help"
  echo ""
  echo "With no arguments - show help"
  exit 1
fi

if [ $# -lt 1 ]; then
    echo "Need 1 parameter! <folder_path>"
    exit 1
fi

# Optional load technical stuff:
source $(dirname $0)/../parallel/util.sh 2> /dev/null

FOLDER=$1
MOD=$2
cd ${FOLDER}

if [ ${COUNT} -gt 0 ]; then
    ALL_COUNT=`expr ${COUNT} - 1`
    OD_COUNT=`expr ${COUNT} - 1`
    YD_COUNT=`expr ${COUNT} - 1`
fi

if [ $(find . -maxdepth 1 -wholename "./*island.bed" | grep -v outlier | wc -l) -gt 0 ]; then
    TOOL="sicer"
fi
if [ $(find . -maxdepth 1 -wholename "./*Peak" | grep -v outlier | wc -l) -gt 0 ]; then
    TOOL="macs2"
fi
if [ $(find . -maxdepth 1 -wholename "./*_peaks.bed" | grep -v outlier | wc -l) -gt 0 ]; then
    TOOL="zinbra"
fi

if [ ${PERCENT} -gt 0 ]; then
    ALL_COUNT=$(echo $(find . -maxdepth 1 \( -wholename "*island.bed" -or -wholename "*Peak" -or \
        -wholename "*_peaks.bed" \) | grep -v outlier | wc -l) " * ${PERCENT} / 100.0 - 1" | bc)
    OD_COUNT=$(echo $(find . -maxdepth 1 \( -wholename "*OD*island.bed" -or -wholename "*OD*Peak" \
        -or -wholename "*OD*_peaks.bed" \) |
        grep -v outlier | wc -l) " * ${PERCENT} / 100.0 - 1" | bc)
    YD_COUNT=$(echo $(find . -maxdepth 1 \( -wholename "*YD*island.bed" -or -wholename "*YD*Peak" \
        -or -wholename "*YD*_peaks.bed" \) |
        grep -v outlier | wc -l) " * ${PERCENT} / 100.0 - 1" | bc)
fi

bedtools multiinter -i $(find . -maxdepth 1 \( -wholename "*_peaks.bed" -or -wholename "*Peak" -or \
    -wholename "*island.bed" \) | grep -v outlier | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/ /g') |\
    grep '\(,.*\)\{'${ALL_COUNT}'\}' | bedtools merge > ${MOD}_${TOOL}_consensus.bed
bedtools multiinter -i $(find . -maxdepth 1 \( -wholename "*OD*island.bed" -or \
    -wholename "*OD*Peak" -or -wholename "*OD*_peaks.bed" \) | grep -v outlier | \
    sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/ /g') | grep '\(,.*\)\{'${OD_COUNT}'\}' | \
    bedtools merge > ${MOD}_${TOOL}_ODS_consensus.bed
bedtools multiinter -i $(find . -maxdepth 1 \( -wholename "*YD*island.bed" -or \
    -wholename "*YD*Peak" -or -wholename "*YD*_peaks.bed" \) | grep -v outlier | \
    sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/ /g') | grep '\(,.*\)\{'${YD_COUNT}'\}' | \
    bedtools merge > ${MOD}_${TOOL}_YDS_consensus.bed
bedtools intersect -v -a ${MOD}_${TOOL}_ODS_consensus.bed \
    -b ${MOD}_${TOOL}_YDS_consensus.bed \
    > ${MOD}_${TOOL}_ODS_without_YDS_consensus.bed
bedtools intersect -v -a ${MOD}_${TOOL}_YDS_consensus.bed \
    -b ${MOD}_${TOOL}_ODS_consensus.bed \
    > ${MOD}_${TOOL}_YDS_without_ODS_consensus.bed
