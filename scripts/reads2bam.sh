#!/usr/bin/env bash
# author: oleg.shpynov@jetbrains.com

# Check tools
which bedtools &>/dev/null || {
    echo "bedtools not found! You can install it using:"
    echo "  conda install -c bioconda bedtools"
    echo "For further details see http://code.google.com/p/bedtools"
    exit 1;
   }


>&2 echo "reads2bam $@"
if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <READS> <CHROM.SIZES>"
    echo "READS: bam, bed, bed.gz"
    exit 1
fi

INPUT=$1
CHROM_SIZES=$2

if [ ! -f "${CHROM_SIZES}" ]; then
  echo "File not found: ${CHROM_SIZES}"
  exit 1
fi

case "$INPUT" in
  *.bed.gz )
    >&2 echo "bed.gz: $INPUT"
    gunzip --keep ${INPUT}
    UNZIPPED=${INPUT%%.gz}
    BAM=${INPUT/.bed.gz/.bam}
    bedtools bedtobam -i ${UNZIPPED} -g ${CHROM_SIZES} > ${BAM}
    rm -f ${UNZIPPED}
    ;;
  *.bed )
    >&2 echo "bed: $INPUT"
    BAM=${INPUT/.bed/.bam}
    bedtools bedtobam -i ${INPUT} -g ${CHROM_SIZES} > ${BAM}
    ;;
  *.bam )
    >&2 echo "bam: $INPUT"
    BAM=${INPUT}
    ;;
  * )
    >&2 echo "UNKNOWN: $INPUT"
    exit 1;
    ;;
esac

echo ${BAM}