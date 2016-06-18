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
#PBS -N subsample_$ID
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o $WORK_DIR/qsub/${NAME}_subsample_${SUBSAMPLE_READS}.log

module load samtools

# This is necessary because qsub default working dir is user home
cd $WORK_DIR

samtools view -h -o ${NAME}.sam ${BAM_FILE}
grep -v '^@' ${NAME}.sam > nohead_${NAME}.sam
sample --sample-size ${SUBSAMPLE_READS} --sample-without-replacement --preserve-order nohead_${NAME}.sam > nohead_${ID}.sam
samtools reheader nohead_${ID}.sam ${BAM_FILE}
samtools view -bS nohead_${ID}.sam -o ${ID}.bam

# Cleanup
rm nohead_$NAME*
rm ${NAME}.sam
ENDINPUT
)
fi