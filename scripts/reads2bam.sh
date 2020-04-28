#!/usr/bin/env bash
# author: oleg.shpynov@jetbrains.com

# Check tools
which bedtools &>/dev/null || {
    echo "ERROR: bedtools not found! You can install it using:"
    echo "  conda install -c bioconda bedtools"
    echo "For further details see http://code.google.com/p/bedtools"
    exit 1;
   }
which samtools &>/dev/null || {
    echo "ERROR: samtools not found! You can install it using:"
    echo "  conda install -c bioconda samtools"
    echo "For further details see http://www.htslib.org/doc/samtools.html"
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

# Check configuration
[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util.sh

export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p "${TMPDIR}"

case "$INPUT" in
  *.bed.gz )
    BAM=${INPUT/.bed.gz/.bam}
    >&2 echo "bed.gz: $INPUT -> ${BAM}"
    if [[ ! -f ${BAM} ]]; then
        UNZIPPED=${INPUT%%.gz}
        gunzip -c ${INPUT} > ${UNZIPPED}
        NS=${BAM}_not_sorted.bam
        bedtools bedtobam -i ${UNZIPPED} -g ${CHROM_SIZES} > ${NS}
        rm -f ${UNZIPPED}
        samtools sort -T ${TMPDIR}/bam.sorted -o ${BAM} ${NS}
        rm ${NS}
    fi
    ;;
  *.bed )
    BAM=${INPUT/.bed/.bam}
    >&2 echo "bed: $INPUT -> ${BAM}"
    if [[ ! -f ${BAM} ]]; then
        NS=${BAM}_not_sorted.bam
        bedtools bedtobam -i ${INPUT} -g ${CHROM_SIZES} > ${NS}
        samtools sort -T ${TMPDIR}/bam.sorted -o ${BAM} ${NS}
        rm ${NS}
    fi
    ;;
  *.bam )
    BAM=${INPUT}
    ;;
  * )
    >&2 echo "UNKNOWN: $INPUT"
    exit 1;
    ;;
esac

# Cleanup
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir

# Result
echo ${BAM}