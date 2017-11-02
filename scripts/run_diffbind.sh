#!/usr/bin/env bash

# Check tools
which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
which macs2 &>/dev/null || { echo "macs2 not found! Install macs2: <https://github.com/taoliu/MACS/wiki/Install-macs2>"; exit 1; }

# Load technical stuff
source $(dirname $0)/../parallel/util.sh
SCRIPT_DIR="$(project_root_dir)"
export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p "${TMPDIR}"

################################################################################
# Configuration start ##########################################################
################################################################################

>&2 echo "run_diffbind: $@"
if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <NAME> <DIFFBIND_CSV>"
    exit 1
fi

NAME=$1
DIFFBIND_CSV=$3

FOLDER=$(pwd)
echo "FOLDER"
echo $FOLDER

PREFIX="$FOLDER/$NAME"
echo "PREFIX"
echo $PREFIX

Q=0.01
echo "Q $Q"
BROAD_CUTOFF=0.1
echo "BROAD_CUTOFF $BROAD_CUTOFF"

READS_Y=$(awk -v FS=',' '{ if ($4 == "Y") print $6 }' $DIFFBIND_CSV | sort  -T ${TMPDIR} --unique | tr '\n' ' ')
echo "READS Y"
echo "$READS_Y"
READS_O=$(awk -v FS=',' '{ if ($4 == "O") print $6 }' $DIFFBIND_CSV | sort -T ${TMPDIR} --unique | tr '\n' ' ')
echo "READS O"
echo "$READS_O"

INPUTS_Y=$(awk -v FS=',' '{ if ($4 == "Y") print $8 }' $DIFFBIND_CSV | sort -T ${TMPDIR} --unique | tr '\n' ' ')
echo "INPUT_READS Y"
echo "$INPUTS_Y"
INPUTS_O=$(awk -v FS=',' '{ if ($4 == "O") print $8 }' $DIFFBIND_CSV | sort -T ${TMPDIR} --unique | tr '\n' ' ')
echo "INPUT_READS O"
echo "$INPUTS_O"

PEAKS_Y=$(awk -v FS=',' '{ if ($4 == "Y") print $9 }' $DIFFBIND_CSV | sed 's#xls#broadPeak#g' | sort -T ${TMPDIR} --unique | tr '\n' ' ')
echo "PEAKS Y"
echo "$PEAKS_Y"
PEAKS_O=$(awk -v FS=',' '{ if ($4 == "O") print $9 }' $DIFFBIND_CSV | sed 's#xls#broadPeak#g' | sort -T ${TMPDIR} --unique | tr '\n' ' ')
echo "PEAKS O"
echo "$PEAKS_O"

################################################################################
# Configuration end ############################################################
################################################################################

DIFFBIND="${PREFIX}_diffbind"
echo
echo "Processing $DIFFBIND"
if [ ! -d $DIFFBIND ]; then
    mkdir -p ${DIFFBIND}
    cp ${DIFFBIND_CSV} ${DIFFBIND}/${NAME}.csv
    cd ${DIFFBIND}

    echo "Processing diffbind"
    run_parallel << SCRIPT
#!/bin/sh
#PBS -N diffbind_${NAME}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=64gb
#PBS -j oe
#PBS -o ${NAME}_diffbind.log
# This is necessary because qsub default working dir is user home
cd ${DIFFBIND}
module load R
Rscript ${SCRIPT_DIR}/R/diffbind.R ${NAME}.csv
SCRIPT
    wait_complete "$QSUB_ID"

    # Filter out old and young donors and sort by Q-Value
    cat ${NAME}_result.csv | awk 'NR > 1 {print $0}' | sed 's#"##g' | tr ',' '\t' |\
        awk -v OFS='\t' '$9 > 0 {print $0}' | sort -T ${TMPDIR} -k10,10g > ${NAME}_cond1.bed
    cat ${NAME}_result.csv | awk 'NR > 1 {print $0}' | sed 's#"##g' | tr ',' '\t' |\
        awk -v OFS='\t' '$9 < 0 {print $0}' | sort -T ${TMPDIR} -k10,10g > ${NAME}_cond2.bed
    # Save ${NAME} results to simple BED3 format
    awk -v OFS='\t' '{print $1,$2,$3}' ${NAME}_cond1.bed > ${NAME}_cond1.bed3
    awk -v OFS='\t' '{print $1,$2,$3}' ${NAME}_cond2.bed > ${NAME}_cond2.bed3
fi

# TMP dir cleanup:
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir
