#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <GENOME> <FOLDER>"
    exit 1
fi
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
    rsync -avzP --exclude="*.txt" rsync://hgdownload.cse.ucsc.edu/goldenPath/${GENOME}/chromosomes/ .
    gunzip *.fa.gz
    chmod a+r *
fi

echo "Check 2bit reference"
if [ ! -f "${GENOME}.2bit" ]; then
    echo "Downloading 2bit reference"
    wget http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME}/bigZips/${GENOME}.2bit
    chmod a+r *
fi

echo "Check chrom.sizes file"
if [ ! -f "${GENOME}.chrom.sizes" ]; then
    echo "Downloading chrom.sizes file"
    wget http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME}/bigZips/${GENOME}.chrom.sizes
    chmod a+r *
fi

echo "Check fa reference"
if [ ! -f "${GENOME}.fa" ]; then
    echo "Convert 2bit to fa reference"
    ~/twoBitToFa/twoBitToFa ${GENOME}.2bit ${GENOME}.fa
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
