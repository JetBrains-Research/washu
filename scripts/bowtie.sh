#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

WORK_DIR=`pwd`
GENOME=$1
FASTQ_FILE=$2
NAME=${FASTQ_FILE%%.f*q} # file name without extension

if [ ! -f "$NAME.bam" ]; then
    echo $(qsub << ENDINPUT
#!/bin/sh
#PBS -N bowtie_${GENOME}_$NAME
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=48gb
#PBS -j oe
#PBS -q dque

# Loading modules
module load bowtie
module load samtools

cd $WORK_DIR
bowtie -p 8 -St -m 1 -v 3 --best --strata $GENOME $FASTQ_FILE $NAME.sam
samtools view -bS -o $NAME.bam $NAME.sam
samtools sort $NAME.bam $NAME.sorted
mv $NAME.sorted.bam $NAME.bam
rm $NAME.sam
ENDINPUT
)
fi