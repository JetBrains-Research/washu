#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com


GENOME=$1
FOLDER=$2

# Load technical stuff
source ~/work/washu/scripts/util.sh

echo "Processing genome $GENOME"
if [ ! -d "$FOLDER" ]; then
    mkdir ${FOLDER}
fi
cd ${FOLDER}

if [ ! -f "chr1.fa" ]; then
    # Download only chromosomes sequences
    rsync -avzP --exclude="chr*_*" --exclude="*.txt" rsync://hgdownload.cse.ucsc.edu/goldenPath/${GENOME}/chromosomes/ .
    gunzip *.fa.gz
    chmod a+r *
fi

echo "Check bowtie indexes ${GENOME}"
if [ ! -f "$GENOME.1.ebwt" ]; then
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N bowtie_indexes_${GENOME}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${FOLDER}/${GENOME}_bowtie_indexes.log

# Load module
module load bowtie

# This is necessary because qsub default working dir is user home
cd ${FOLDER}
bowtie-build $(find . -type f -name "*.fa" -printf '%P\n' | paste -sd "," -) ${GENOME}
ENDINPUT
)
    wait_complete ${QSUB_ID}
    check_logs
fi