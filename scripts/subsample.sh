#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

WORK_DIR=`pwd`
BAM_FILE=$1
SUBSAMPLE_READS=$2
NAME=${BAM_FILE%%.bam} # file name without extension
ID=${NAME}_${SUBSAMPLE_READS}

if [ ! -f "${ID}.bam" ]; then
    echo $(qsub << ENDINPUT
#!/bin/sh
#PBS -N subsample_${ID}
#PBS -l nodes=1:ppn=8,walltime=4:00:00,vmem=8gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_subsample_${SUBSAMPLE_READS}.log

module load samtools

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

samtools view -H ${BAM_FILE} > head_${NAME}.sam
samtools view ${BAM_FILE} > nohead_${NAME}.sam
shuf -n ${SUBSAMPLE_READS} nohead_${NAME}.sam > nohead_${ID}.sam
cat head_${NAME}.sam nohead_${ID}.sam > ${ID}.sam
samtools view -bS ${ID}.sam -o not_sorted_${ID}.bam
samtools sort not_sorted_${ID}.bam -o ${ID}.bam

# Cleanup
rm *head*_${NAME}*
rm ${ID}.sam
rm not_sorted_${ID}.bam
ENDINPUT
)
fi