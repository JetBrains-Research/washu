#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

WORK_DIR=`pwd`
FASTQ_FILE=$1
NAME=${FASTQ_FILE%%.f*q} # file name without extension

# Fastq produce 2 files: _fastqc.html and _fastq.zip
if [ ! -f "${NAME}_fastqc.html" ]; then
    echo $(qsub << ENDINPUT
#!/bin/sh
#PBS -N fastq_$NAME
#PBS -l nodes=1:ppn=8,walltime=2:00:00,vmem=6gb
#PBS -j oe
#PBS -q dque

# Loading modules
module load fastqc

cd $WORK_DIR
fastqc $NAME
ENDINPUT
)
fi