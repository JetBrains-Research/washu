#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Load technical stuff
source ~/work/washu/scripts/util.sh

WORK_DIR=$1
GENOME=$2
FOLDER=$3
Q=$4

echo "Batch zinbra: ${WORK_DIR} ${GENOME} ${FOLDER} ${Q}"
cd ${WORK_DIR}

TASKS=""
for FILE in $(find . -type f -name '*.bam' -printf '%P\n')
do :
    NAME=${FILE%%.bam} # file name without extension
    ID=${NAME}_${GENOME}_${Q}

    # Submit task
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N zinbra_${ID}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_zinbra_${GENOME}.log

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

module load java
export _JAVA_OPTIONS="-Xms512m -Xmx14g"
java -jar /home/oshpynov/zinbra/zinbra-0.2.4.jar analyze --input ${FILE} --reference ${FOLDER}/${GENOME}.2bit --fdr ${Q} --bed ${ID}_peaks.bed

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

echo "Done. Batch zinbra: ${WORK_DIR} ${GENOME} ${FOLDER} ${Q}"