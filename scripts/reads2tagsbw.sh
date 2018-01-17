#!/usr/bin/env bash
# author: oleg.shpynov@jetbrains.com

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
source ${WASHU_ROOT}/parallel/util/util.sh/null


# Convert reads to BAM is required
BAM=$(bash "${WASHU_ROOT}/scripts/reads2bam.sh" ${INPUT} ${CHROM_SIZES})
# Covert bam to bdg
if [[ ! -z  ${RESULT} ]]; then
    BDG=${RESULT/.bw/.bdg}
else
    BDG=${BAM/.bam/_tags.bdg}
fi

export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p "${TMPDIR}"

NAME=${BAM%%.bam}; NAME=${NAME##*/} # file name without extension

bash ${WASHU_ROOT}/scripts/bam2tags.sh ${BAM} ${INSERT_SIZE} > ${TMPDIR}/${NAME}.tags
bash ${WASHU_ROOT}/scripts/tags2bdg.sh  ${TMPDIR}/${NAME}.tags > ${BDG}
bash ${WASHU_ROOT}/scripts/bdg2bw.sh ${BDG} ${CHROM_SIZES}
# Cleanup
rm ${BDG}

# TMP dir cleanup:
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir