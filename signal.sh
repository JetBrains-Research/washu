#!/usr/bin/env bash

# See MACS2 docs for signal track: https://github.com/taoliu/MACS/wiki/Build-Signal-Track
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
bash ~/work/washu/bdg2bw.sh ${ID}_FE.bdg ${CHROM_SIZES}
