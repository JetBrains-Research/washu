#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Load technical stuff
source $(dirname $0)/../parallel/util.sh

>&2 echo "Batch subsample $@"
if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <WORK_DIR> <READS>"
    exit 1
fi
WORK_DIR=$1
READS=$2

cd ${WORK_DIR}

TASKS=""
for FILE in $(find . -name '*.bam' | sed 's#\./##g')
do :
    NAME=${FILE%%.bam} # file name without extension
    ID=${NAME}_${READS}mln

    # Submit task
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N subsample_${ID}
#PBS -l nodes=1:ppn=8,walltime=4:00:00,vmem=64gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${ID}.log

module load samtools

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

samtools view -H ${FILE} > head_${NAME}.sam
samtools view ${FILE} > nohead_${NAME}.sam
shuf -n ${READS}000000 nohead_${NAME}.sam > nohead_${ID}.sam
cat head_${NAME}.sam nohead_${ID}.sam > ${ID}.sam
samtools view -bS ${ID}.sam -o not_sorted_${ID}.bam
samtools sort not_sorted_${ID}.bam -o ${ID}.bam

# Cleanup
rm *head*_${NAME}*
rm ${ID}.sam
rm not_sorted_${ID}.bam
ENDINPUT
)
    echo "FILE: ${FILE}; TASK: ${QSUB_ID}"
    TASKS="$TASKS $QSUB_ID"
done
wait_complete ${TASKS}
check_logs

>&2 echo "Done. Batch subsample $@"