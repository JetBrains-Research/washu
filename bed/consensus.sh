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
TOOL=""
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
    echo "Need 1 parameter! <files>"
    exit 1
fi

# Optional load technical stuff:
source $(dirname $0)/../parallel/util/util.sh 2> /dev/null

FILES=$@

if [ ${COUNT} -gt 0 ]; then
    CONS_COUNT=`expr ${COUNT} - 1`
fi

if [[ ${FILES} == *"island.bed"* ]]; then
  TOOL=${TOOL}"_sicer"
fi
if [[ ${FILES} == *"Peak"* ]]; then
  TOOL=${TOOL}"_macs2"
fi
if [[ ${FILES} == *"_peaks.bed"* ]]; then
  TOOL=${TOOL}"_zinbra"
fi

if [ ${PERCENT} -gt 0 ]; then
    CONS_COUNT=$(echo "(" $(echo "${FILES}" | awk -F" " '{print NF-1}') \
        " + 1) * ${PERCENT} / 100.0 - 1" | bc)
fi

bedtools multiinter -i ${FILES} | grep '\(,.*\)\{'${CONS_COUNT}'\}' | \
    sort -k1,1 -k2,2n | bedtools merge
