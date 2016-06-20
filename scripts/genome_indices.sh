#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com


GENOME=$1
WORK_DIR=$2

# Load technical stuff
source ~/work/washu/scripts/util.sh

echo "Processing genome $GENOME"
if [ ! -d "$WORK_DIR/$GENOME" ]; then
    mkdir $WORK_DIR/$GENOME
fi
if [ ! -f "$WORK_DIR/$GENOME/chr1.fa" ]; then
    cd $WORK_DIR/$GENOME
    # Download only chromosomes sequences
    rsync -avzP --exclude="chr*_*" --exclude="*.txt" rsync://hgdownload.cse.ucsc.edu/goldenPath/$GENOME/chromosomes/ .
    gunzip *.fa.gz
    chmod a+r *
    cd $WORK_DIR
fi

echo "Check bowtie indexes $GENOME"
if [ ! -f "$WORK_DIR/$GENOME/$GENOME.1.ebwt" ]; then
    cd $WORK_DIR/$GENOME
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N bowtie_indexes_${GENOME}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o $WORK_DIR/qsub/bowtie_indexes_${GENOME}.log

# Load module
module load bowtie

# This is necessary because qsub default working dir is user home
cd $WORK_DIR/$GENOME
bowtie-build $(find $WORK_DIR/$GENOME -type f -name "*.fa" -printf '%P\n' | paste -sd "," -) $GENOME
ENDINPUT
)
    cd $WORK_DIR
    wait_complete $QSUB_ID
    check_logs
fi