#!/usr/bin/env bash
# Script to create report for signals
# Output example:
# H3K27me3  150  cpg_minavcov10_complex_4outliers.narrow.adjusted.regions.filtered  rawz                              2
# H3K27ac   150  cpg_minavcov10_complex_4outliers.narrow.adjusted.regions.filtered  diffbind_tmm_minus_full           3
# author oleg.shpynov@jetbrains.com

>&2 echo "signals_report.sh $@"
if [ $# -lt 1 ]; then
    echo "Need 1 parameter! <WORK_DIR>"
    exit 1
fi

WORK_DIR=$1

cd $WORK_DIR
T=$'\t';
echo "modification${T}fragment${T}file${T}normalization${T}e"
for F in $(find . -name "*_pca_fit_error.csv"); do
    IFS='/' read -ra CHUNKS <<< "$F"
    MODIFICATION=${CHUNKS[1]}
    FRAGMENT=${CHUNKS[2]}
    FOLDER=${CHUNKS[3]}
    FILE=${CHUNKS[4]}
    NORMALIZATION=$(echo "$FILE" | sed 's#_pca_fit_error.csv##g' | sed "s#${FOLDER}_##g" | sed "s#.*\]_##g")
    ERROR=$(cat $F)
    echo "$MODIFICATION$T$FRAGMENT$T$FOLDER$T$NORMALIZATION$T$ERROR"
done