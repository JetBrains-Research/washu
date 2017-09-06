#!/usr/bin/env bash

# Original code: http://crazyhottommy.blogspot.ru/2014/10/convert-bam-file-to-bigwig-file-and.html
# author: oleg.shpynov@jetbrains.com

# Check tools
which bedtools &>/dev/null || {
    echo "bedtools not found! You can install it using:"
    echo "  conda install -c bioconda bedtools"
    echo "For further details see http://code.google.com/p/bedtools"
    exit 1;
   }

# Load technical stuff
source $(dirname $0)/../parallel/util.sh
SCRIPT_DIR="$(project_root_dir)"

>&2 echo "bam2bw $@"
if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <BAM> <chrom.sizes>"
    exit 1
fi

BAM=$1
CHROM_SIZES=$2

if [ ! -f "${CHROM_SIZES}" ]; then
  echo "File not found: ${CHROM_SIZES}"
  exit 1
fi

NAME=${BAM%%.bam}
bedtools genomecov -ibam $BAM -bg -g ${CHROM_SIZES} > ${NAME}.bdg
bash "${SCRIPT_DIR}/scripts/bdg2bw.sh" ${NAME}.bdg ${CHROM_SIZES}