#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

WORK_DIR=`pwd`
SRA_FILE=$1
NAME=${SRA_FILE%%.sra} # file name without sra extension

if [ ! -f "${NAME}.fq" ]; then
    echo $(qsub << ENDINPUT
#!/bin/sh
#PBS -N sra2fastq_$NAME
#PBS -l nodes=1:ppn=1,walltime=2:00:00,vmem=8gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_sra2fastq.log

# Loading sratoolkit module
module load sratoolkit

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}
fastq-dump --split-3 --outdir ${WORK_DIR} ${SRA_FILE}
ENDINPUT
)
fi