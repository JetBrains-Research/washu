#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <GENOME> <FOLDER>"
    exit 1
fi
GENOME=$1
FOLDER=$2

# Load technical stuff, not available in qsub emulation
if [ -f "$(dirname $0)/util.sh" ]; then
    source "$(dirname $0)/util.sh"
fi

echo "Processing genome $GENOME"
if [ ! -d "$FOLDER" ]; then
    mkdir -p ${FOLDER}
fi
cd ${FOLDER}

if [ ! -f "chr1.fa" ]; then
    # Download only chromosomes sequences
    rsync -avzP --exclude="*.txt" rsync://hgdownload.cse.ucsc.edu/goldenPath/${GENOME}/chromosomes/ .
    gunzip *.fa.gz
    chmod a+r *
fi

echo "Check chrom.sizes file"
if [ ! -f "${GENOME}.chrom.sizes" ]; then
    echo "Downloading chrom.sizes file"
    wget http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME}/bigZips/${GENOME}.chrom.sizes
    chmod a+r *
fi
