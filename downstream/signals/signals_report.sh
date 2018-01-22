#!/usr/bin/env bash
# Script to create report for signals
# author oleg.shpynov@jetbrains.com

>&2 echo "signals_report.sh $@"
if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <WORK_DIR>"
    exit 1
fi

WORK_DIR=$1

cd $WORK_DIR
T=$'\t';
echo "modification${T}file${T}normalization${T}e"
for F in $(find . -name "*_fit_error.csv"); do
    IFS='/' read -ra CHUNKS <<< "$F"
    MODIFICATION=${CHUNKS[1]}
    FOLDER=${CHUNKS[2]}
    FILE=${CHUNKS[3]}
    NORMALIZATION=$(echo "$FILE" | sed 's#_fit_error.csv##g' | sed "s#${FOLDER}_##g")
    ERROR=$(cat $F)
    echo "$MODIFICATION$T$FOLDER$T$NORMALIZATION$T$ERROR"
done