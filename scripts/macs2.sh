#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Load technical stuff
source ~/work/washu/scripts/util.sh

WORK_DIR=$1
GENOME=$2
Q=$3

echo "Batch macs2: ${WORK_DIR} ${GENOME} ${Q}"

# Convert Genome build to macs2 species
[[ ${GENOME} =~ ^hg[0-9]+$ ]] && SPECIES="hs"
[[ ${GENOME} =~ ^mm[0-9]+$ ]] && SPECIES="mm"
[[ -z "$SPECIES" ]] && echo "Unknown species for macs: $GENOME" && exit 1

cd ${WORK_DIR}

TASKS=""
for FILE in $(find . -type f -name '*.bam' -printf '%P\n')
do :
    NAME=${FILE%%.bam} # file name without extension
    ID=${NAME}_${GENOME}_${Q}

    # Submit task
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N macs2_${ID}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_macs2_${GENOME}.log

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}
/home/oshpynov/miniconda2/bin/macs2 callpeak -t ${FILE} -f BAM -g ${SPECIES} -n ${ID} -B -q ${Q}

# Cleanup
mv ${ID}_peaks.narrowPeak do_not_remove_${ID}_peaks.bed
rm ${ID}*
mv do_not_remove_${ID}_peaks.bed ${ID}_peaks.bed
ENDINPUT
)
    echo "$FILE: $QSUB_ID"
    TASKS="$TASKS $QSUB_ID"
done
wait_complete ${TASKS}
check_logs

echo "Batch macs2: ${WORK_DIR} ${GENOME} ${Q}"