#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Load technical stuff
source $(dirname $0)/../parallel/util.sh

>&2 echo "index-genome $@"
if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <GENOME> <FOLDER>"
    exit 1
fi
GENOME=$1
FOLDER=$2


if [ ! -d "$FOLDER" ]; then
    mkdir -p ${FOLDER}
fi
cd ${FOLDER}

if [ ! -f "chr22.fa" ]; then
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

>&2 echo "Done. index-genome $@"
