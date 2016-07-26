#!/usr/bin/env bash

# Original code: https://github.com/ENCODE-DCC/chip-seq-pipeline/blob/master/dnanexus/macs2/src/macs2.py#L143
# author oleg.shpynov@jetbrains.com

which macs2 &>/dev/null || { echo "macs2 not found!"; exit 1; }
which bedtools &>/dev/null || { echo "bedtools not found!"; exit 1; }

if [ $# -lt 4 ]; then
    echo "Need 4 parameters! <ID> <TREAT> <CONTROL> <chrom.sizes>"
    exit 1
fi

ID=$1
TREAT=$2
CONTROL=$3
CHROM_SIZES=$4

if [ ! -f "${CHROM_SIZES}" ]; then
  echo "File not found: ${CHROM_SIZES}"
  exit 1
fi

echo "Create fold enrichment signal track"
macs2 bdgcmp -t ${TREAT} -c ${CONTROL} -o ${ID}_FE.bdg -m FE

# Remove coordinates outside chromosome sizes (stupid MACS2 bug)
slopBed -i ${ID}_FE.bdg -g ${ID} -b 0 | bedClip stdin ${ID} ${ID}_signal.bedgraph
# To BigWig
bedGraphToBigWig ${ID}_signal.bedgraph ${CHROM_SIZES} ${ID}_signal.bw