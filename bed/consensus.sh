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

if [ ${COUNT} -gt 0 ]; then
    ALL_COUNT=${COUNT}
    OD_COUNT=${COUNT}
    YD_COUNT=${COUNT}
fi

# comments about percent and count
for WORK_DIR in $(find $FOLDER -type d -name "H*"); do :
    WORK_DIR_NAME=${WORK_DIR##*/}
    echo  "!${WORK_DIR_NAME}"

    cd ${WORK_DIR}

    if [ $(find . -wholename "./bed/*island.bed" | grep -v outlier | wc -l) -gt 0 ]; then
        TOOL="sicer"
    fi
    if [ $(find . -wholename "./bed/*Peak" | grep -v outlier | wc -l) -gt 0 ]; then
        TOOL="macs2"
    fi
    if [ $(find . -wholename "./bed/*_peaks.bed" | grep -v outlier | wc -l) -gt 0 ]; then
        TOOL="zinbra"
    fi

    if [ ${PERCENT} -gt 0 ]; then
        ALL_COUNT=$(echo $(find . \( -wholename "./bed/*island.bed" -or -wholename "./bed/*Peak" -or
            -wholename "./bed/*_peaks.bed" \) |
            grep -v outlier | wc -l) * \(${PERCENT} / 100.0\) - 1 | bc)
        OD_COUNT=$(echo $(find . \( -wholename "./bed/*OD*island.bed" -or
            -wholename "./bed/*OD*Peak" -or -wholename "./bed/*OD*_peaks.bed" \) |
            grep -v outlier | wc -l) * \(${PERCENT} / 100.0\) - 1 | bc)
        YD_COUNT=$(echo $(find . \( -wholename "./bed/*YD*island.bed" -or
            -wholename "./bed/*YD*Peak" -or -wholename "./bed/*YD*_peaks.bed" \) |
            grep -v outlier | wc -l) * \(${PERCENT} / 100.0\) - 1 | bc)
    fi

    find . \( -wholename "./bed/*island.bed" -or -wholename "./bed/*Peak" -or
        -wholename "./bed/*_peaks.bed" \) | grep -v outlier |
        sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/ /g' |
        xargs -0 bash /mnt/stripe/washu/bed/union.sh | grep '\(|.*\)\{'${ALL_COUNT}'\}' |
        awk -v OFS='\t' '{print $1,$2,$3}' > ${WORK_DIR_NAME}_${TOOL}_median_consensus.bed
    find . \( -wholename "./bed/*OD*island.bed" -or -wholename "./bed/*OD*Peak" -or
        -wholename "./bed/*OD*_peaks.bed" \) | grep -v outlier |
        sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/ /g' | xargs -0 bash /mnt/stripe/washu/bed/union.sh |
        grep '\(|.*\)\{'${OD_COUNT}'\}' | awk -v OFS='\t' '{print $1,$2,$3}'
        > ${WORK_DIR_NAME}_${TOOL}_ODS_median_consensus.bed
    find . \( -wholename "./bed/*YD*island.bed" -or -wholename "./bed/*YD*Peak" -or
        -wholename "./bed/*YD*_peaks.bed" \) | grep -v outlier |
        sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/ /g' | xargs -0 bash /mnt/stripe/washu/bed/union.sh |
        grep '\(|.*\)\{'${YD_COUNT}'\}' | awk -v OFS='\t' '{print $1,$2,$3}'
        > ${WORK_DIR_NAME}_${TOOL}_YDS_median_consensus.bed
    bedtools intersect -v -a ${WORK_DIR_NAME}_${TOOL}_ODS_median_consensus.bed -b
        ${WORK_DIR_NAME}_${TOOL}_YDS_median_consensus.bed
        > ${WORK_DIR_NAME}_${TOOL}_ODS_without_YDS_median_consensus.bed
    bedtools intersect -v -a ${WORK_DIR_NAME}_${TOOL}_YDS_median_consensus.bed -b
        ${WORK_DIR_NAME}_${TOOL}_ODS_median_consensus.bed
        > ${WORK_DIR_NAME}_${TOOL}_YDS_without_ODS_median_consensus.bed
    #echo $(pwd)
done