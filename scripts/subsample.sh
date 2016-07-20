#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Load technical stuff
source ~/work/washu/scripts/util.sh

WORK_DIR=$1
READS=$2

echo "Batch subsampling: ${WORK_DIR} ${READS}"
cd ${WORK_DIR}

TASKS=""
for FILE in $(find . -type f -name '*.bam' -printf '%P\n')
do :
    NAME=${FILE%%.bam} # file name without extension
    ID=${NAME}_${READS}

    # Submit task
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N subsample_${ID}
#PBS -l nodes=1:ppn=8,walltime=4:00:00,vmem=8gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_subsample_${READS}.log

module load samtools

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

samtools view -H ${FILE} > head_${NAME}.sam
samtools view ${FILE} > nohead_${NAME}.sam
shuf -n ${READS} nohead_${NAME}.sam > nohead_${ID}.sam
cat head_${NAME}.sam nohead_${ID}.sam > ${ID}.sam
samtools view -bS ${ID}.sam -o not_sorted_${ID}.bam
samtools sort not_sorted_${ID}.bam -o ${ID}.bam

# Cleanup
rm *head*_${NAME}*
rm ${ID}.sam
rm not_sorted_${ID}.bam
ENDINPUT
)
    echo "FILE: ${FILE}; JOB: ${QSUB_ID}"
    TASKS="$TASKS $QSUB_ID"
done
wait_complete ${TASKS}
check_logs

echo "Done. Batch subsampling: ${WORK_DIR} ${READS}"