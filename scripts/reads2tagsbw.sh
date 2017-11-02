#!/usr/bin/env bash
# author: oleg.shpynov@jetbrains.com

# Check tools
which bedtools &>/dev/null || {
    echo "bedtools not found! You can install it using:"
    echo "  conda install -c bioconda bedtools"
    echo "For further details see http://code.google.com/p/bedtools"
    exit 1;
   }

>&2 echo "reads2tagsbw $@"
if [ $# -lt 3 ]; then
    echo "Need at least 3 parameters! <READS> <INSERT_SIZE> <CHROM_SIZES> [<OUTPUT.bw>]"
    echo "READS: bam, bed, bed.gz"
    exit 1
fi

INPUT=$1
INSERT_SIZE=$2
CHROM_SIZES=$3
RESULT=$4

if [ ! -f "${CHROM_SIZES}" ]; then
  echo "File not found: ${CHROM_SIZES}"
  exit 1
fi

SHIFT=$(($INSERT_SIZE / 2))

# Optional load technical stuff:
source $(dirname $0)/../parallel/util.sh 2> /dev/null
SCRIPT_DIR="$(project_root_dir)"


# Convert reads to BAM is required
BAM=$(bash "${SCRIPT_DIR}/scripts/reads2bam.sh" ${INPUT} ${CHROM_SIZES})
NAME=${BAM%%.bam}; NAME=${NAME##*/} # file name without extension
# Covert bam to bdg
if [[ ! -z  ${RESULT} ]]; then
    BDG=${RESULT/.bw/.bdg}
else
    BDG=${BAM/.bam/_tags.bdg}
fi

export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p "${TMPDIR}"

bash ${SCRIPT_DIR}/scripts/bam2tags.sh ${BAM} ${INSERT_SIZE} > ${TMPDIR}/${NAME}.tags
bedtools merge -i ${TMPDIR}/${NAME}.tags -c 1 -o count > ${BDG}
bash ${SCRIPT_DIR}/scripts/bdg2bw.sh ${BDG} ${CHROM_SIZES}
# Cleanup
rm ${BDG}

# TMP dir cleanup:
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir