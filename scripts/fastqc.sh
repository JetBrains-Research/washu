#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

WORK_DIR=`pwd`
FASTQ_FILE=$1
NAME=${FASTQ_FILE%%.f*q} # file name without extension

# Fastq produce 2 files: _fastqc.html and _fastq.zip
if [ ! -f "${NAME}_fastqc.html" ]; then
    echo $(qsub << ENDINPUT
#!/bin/sh
#PBS -N fastqc_${NAME}
#PBS -l nodes=1:ppn=8,walltime=2:00:00,vmem=8gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_fastqc.log

# Loading modules
module load fastqc

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}
fastqc ${FASTQ_FILE}
ENDINPUT
)
fi