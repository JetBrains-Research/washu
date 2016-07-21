#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Load technical stuff
source ~/work/washu/scripts/util.sh

WORK_DIR=$1
GENOME=$2
INDEXES=$3

echo "Batch Bowtie: ${WORK_DIR} ${GENOME} ${INDEXES}"
cd ${WORK_DIR}

TASKS=""
for FILE in $(find . -type f -name '*.f*q' -printf '%P\n')
do :
    NAME=${FILE%%.f*q} # file name without extension
    ID=${NAME}_${GENOME}

    # Submit task
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N bowtie_${GENOME}_${NAME}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_bowtie_${GENOME}.log

# Loading modules
module load bowtie
module load samtools

export BOWTIE_INDEXES=${INDEXES}
# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}
bowtie -p 8 -St -m 1 -v 3 --best --strata ${GENOME} ${FILE} ${ID}.sam
samtools view -bS ${ID}.sam -o ${ID}_not_sorted.bam
samtools sort ${ID}_not_sorted.bam -o ${ID}.bam

# Cleanup
rm ${ID}.sam ${ID}_not_sorted.bam

ENDINPUT
)
    echo "FILE: ${FILE}; JOB: ${QSUB_ID}"
    TASKS="$TASKS $QSUB_ID"
done
wait_complete ${TASKS}
check_logs

echo "Done. Batch Bowtie: ${WORK_DIR} ${GENOME} ${INDEXES}"