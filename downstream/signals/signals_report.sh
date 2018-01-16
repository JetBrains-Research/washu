#!/usr/bin/env bash
# Script to create report for signals
# author oleg.shpynov@jetbrains.com

>&2 echo "signals_report.sh $@"
if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <WORK_DIR> <REPORT_TSV>"
    exit 1
fi

WORK_DIR=$1
REPORT_TSV=$2

# Load technical stuff
source $(dirname $0)/../../parallel/util/util.sh
export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p $TMPDIR
REPORT_TSV_TMP=$(mktemp $TMPDIR/reportXXXXXXX.tsv)

T=$'\t';
cd $WORK_DIR
for F in $(find . -name "*_fit_error.csv"); do
    N=$(echo $F | sed 's#\./##g');
    M=${N%%/*};
    R=${N##*/};
    echo $M; echo $R;
    L=$(cat $F | tr ',' '\t');
    echo "$M$T$R$T$L" >> ${REPORT_TSV_TMP};
done
echo "modification${T}file${T}e${T}e_scaled${T}e_log${T}e_scaled_log${T}e_min" > ${REPORT_TSV}
cat ${REPORT_TSV_TMP} | awk -v OFS='\t' '{min=$3; for(j=4;j<=6;j++){min=($j<min)?$j:min}; print($1,$2,$3,$4,$5,$6,min)}' >> ${REPORT_TSV}

# Cleanup
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir

echo "REPORT for DIR ${WORK_DIR} created: ${REPORT_TSV}"