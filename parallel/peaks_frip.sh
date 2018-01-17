#!/usr/bin/env bash
# This script is used to filter peaks by given FDR from given MACS2 peaks folder
# author Oleg Shpynov (oleg.shpynov@jetbrains.com)
# author Roman Chernyatchik (roman.chernyatchik@jetbrains.com)

which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }

# Load technical stuff
source ${WASHU_ROOT}/parallel/util/util.sh

>&2 echo "peaks_frip.sh: $@"

if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <COMMA_SEPARATED_PEAKS_FOLDERS> <COMMA_SEPARATED_READS_FOLDERS>"
    exit 1
fi

PEAKS_FOLDERS_ARG=$1
READS_FOLDERS_ARG=$2

echo "Compute FRIPs for READS_FOLDER: $READS_FOLDERS_ARG; PEAKS_FOLDER: $PEAKS_FOLDERS_ARG"

PEAKS_FOLDERS=($(echo ${PEAKS_FOLDERS_ARG} | awk -F, '{ for (i=1;i<=NF;i++)print $i }'))
READS_FOLDERS=($(echo ${READS_FOLDERS_ARG} | awk -F, '{ for (i=1;i<=NF;i++)print $i }'))

if [ ${#PEAKS_FOLDERS[*]} -ne ${#READS_FOLDERS[*]} ]; then
  echo "Peaks and reads folders number should be same: ${#PEAKS_FOLDERS[*]} != ${#READS_FOLDERS[*]}"
  exit 1
fi

TASKS=()
for i in ${!PEAKS_FOLDERS[*]}; do
    PEAKS_FOLDER=${PEAKS_FOLDERS[$i]}
    READS_FOLDER=${READS_FOLDERS[$i]}

    if [[ ! -d ${PEAKS_FOLDER} ]]; then
        echo "Missing folder ${PEAKS_FOLDER}"
        exit 1
    fi
    if [[ ! -d ${READS_FOLDER} ]]; then
        echo "Missing folder ${READS_FOLDER}"
        exit 1
    fi

    cd ${PEAKS_FOLDER}
    PEAKS_FOLDER_NAME=${PEAKS_FOLDER##*/}
    for F in $(ls *.*Peak | grep -v gapped); do
        NAME=${F%%_broad*} # filter if broad peaks
        NAME=${NAME%%_q0.*}   # filter suffix if narrow peaks
        BAM=${READS_FOLDER}/${NAME}*.bam

        # Submit task
        run_parallel << SCRIPT
#!/bin/sh
#PBS -N frip_${PEAKS_FOLDER_NAME}_${F}
#PBS -l nodes=1:ppn=1,walltime=2:00:00,vmem=4gb
#PBS -j oe
#PBS -o ${PEAKS_FOLDER}/${F}_frip.log

cd ${PEAKS_FOLDER}
bash ${WASHU_ROOT}/scripts/rip.sh ${BAM} ${F}

SCRIPT
        echo "FILE: ${F}; TASK: ${QSUB_ID}"
        TASKS+=("$QSUB_ID")
    done
done

wait_complete ${TASKS[@]}
check_logs

for PEAKS_FOLDER in ${PEAKS_FOLDERS[*]}; do
    python ${WASHU_ROOT}/parallel/util/peaks_logs.py ${PEAKS_FOLDER}
done