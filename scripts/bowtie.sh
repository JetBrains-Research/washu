#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

WORK_DIR=`pwd`
GENOME=$1
FASTQ_FILE=$2
INDICES=$3
NAME=${FASTQ_FILE%%.f*q} # file name without extension
ID=${NAME}_${GENOME}

if [ ! -f ${ID}.bam ]; then
    echo $(qsub << ENDINPUT
#!/bin/sh
#PBS -N bowtie_${GENOME}_${NAME}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_bowtie_${GENOME}.log

# Loading modules
module load bowtie
module load samtools

export BOWTIE_INDEXES=${INDICES}
# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}
bowtie -p 8 -St -m 1 -v 3 --best --strata ${GENOME} ${FASTQ_FILE} ${ID}.sam
samtools view -bS ${ID}.sam -o ${ID}_not_sorted.bam
samtools sort ${ID}_not_sorted.bam -o ${ID}.bam
# Cleanup
rm ${ID}.sam ${ID}_not_sorted.bam
ENDINPUT
)
fi