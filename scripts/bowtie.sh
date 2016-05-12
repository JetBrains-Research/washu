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
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o $WORK_DIR/qsub/bowtie_${GENOME}_$NAME.log

# Loading modules
module load bowtie
module load samtools

export BOWTIE_INDEXES="$WORK_DIR/$GENOME"
# This is necessary because qsub default working dir is user home
cd $WORK_DIR
bowtie -p 8 -St -m 1 -v 3 --best --strata $GENOME $FASTQ_FILE $NAME.sam
samtools view -bS -o ${NAME}_not_sorted.bam $NAME.sam
samtools sort ${NAME}_not_sorted.bam -o $NAME.bam
# Remove intermediate files
rm $NAME.sam ${NAME}_not_sorted.bam
ENDINPUT
)
fi