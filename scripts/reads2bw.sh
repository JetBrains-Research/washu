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

# Check configuration
[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util.sh

>&2 echo "bam2bw $@"
if [[ $# -lt 2 ]]; then
    echo "Need at least 2 parameters! <READS> <chrom.sizes> [<genes.gtf>]"
    echo "READS: bam, bed, bed.gz"
    exit 1
fi

INPUT=$1
CHROM_SIZES=$2
GTF_FILE=$3

if [[ ! -f "${CHROM_SIZES}" ]]; then
  echo "File not found: ${CHROM_SIZES}"
  exit 1
fi

# Convert reads to BAM is required
BAM=$(bash "${WASHU_ROOT}/scripts/reads2bam.sh" ${INPUT} ${CHROM_SIZES})
# Get correct bedGraph file name
BDG_FILE=$(echo ${BAM} | sed -r 's#\.(bam|bed|bed\.gz)$#.bdg#g')

echo "Converting ${BAM} -> ${BDG_FILE}"
if [[ -f  ${GTF_FILE} ]]; then
    echo "GTF file found, processing paired-end exome alignment"
    bedtools genomecov -ibam ${BAM} -split -bg -g ${CHROM_SIZES} > ${BDG_FILE}
else
    bedtools genomecov -ibam ${BAM} -bg -g ${CHROM_SIZES} > ${BDG_FILE}
fi

# Covert bedGraph to BigWig
bash "${WASHU_ROOT}/scripts/bdg2bw.sh" ${BDG_FILE} ${CHROM_SIZES}

# Cleanup
rm ${BDG_FILE}