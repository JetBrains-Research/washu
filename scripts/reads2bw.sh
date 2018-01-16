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
source $(dirname $0)/../parallel/util/util.sh
SCRIPT_DIR="$(project_root_dir)"

>&2 echo "bam2bw $@"
if [ $# -lt 2 ]; then
    echo "Need at least 2 parameters! <READS> <chrom.sizes> [<OUTPUT.bw>]"
    echo "READS: bam, bed, bed.gz"
    exit 1
fi

INPUT=$1
CHROM_SIZES=$2
RESULT=$3

if [ ! -f "${CHROM_SIZES}" ]; then
  echo "File not found: ${CHROM_SIZES}"
  exit 1
fi

# Convert reads to BAM is required
BAM=$(bash "${SCRIPT_DIR}/scripts/reads2bam.sh" ${INPUT} ${CHROM_SIZES})

# Covert bam to bdg
if [[ ! -z  ${RESULT} ]]; then
    BDG_FILE=${RESULT/.bw/.bdg}
else
    BDG_FILE=${BAM/.bam/.bdg}
fi

bedtools genomecov -ibam $BAM -bg -g ${CHROM_SIZES} > ${BDG_FILE}
bash "${SCRIPT_DIR}/scripts/bdg2bw.sh" ${BDG_FILE} ${CHROM_SIZES}
# Cleanup
rm ${BDG_FILE}