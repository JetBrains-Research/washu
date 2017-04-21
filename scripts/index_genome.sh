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

# NOTE[shpynov]: be careful, since full fa reference can be used for bowtie indexes by mistake,
# which leads to huge multimapped fraction (>99%) in real cases.
#
# Easiest way to obtain full fa file is download 2bit and convert
#echo "Check fa reference"
#if [ ! -f "${GENOME}.fa" ]; then
#    echo "Downloading 2bit reference"
#    wget http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME}/bigZips/${GENOME}.2bit
#    chmod a+r *
#    echo "Convert 2bit to fa reference"
#    ~/twoBitToFa/twoBitToFa ${GENOME}.2bit ${GENOME}.fa
#    chmod a+r *
#fi