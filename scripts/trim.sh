#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

WORK_DIR=`pwd`
FASTQ_FILE=$1
TRIM_START=$2
NAME=${FASTQ_FILE%%.f*q} # file name without extension
ID=${NAME}_${TRIM_START}
if [ ! -f "${ID}.fq" ]; then
    echo $(qsub << ENDINPUT
#!/bin/sh
#PBS -N bowtie_${GENOME}_${NAME}
#PBS -l nodes=1:ppn=8,walltime=4:00:00,vmem=4gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_trim_${TRIM_START}.log

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

/home/oshpynov/seqtk/seqtk trimfq -b 5 ${FASTQ_FILE} > ${ID}.fq
ENDINPUT
)
fi