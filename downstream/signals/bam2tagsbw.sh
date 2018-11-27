#!/usr/bin/env bash
# author: oleg.shpynov@jetbrains.com

>&2 echo "bam2tagsbw $@"
if [ $# -lt 4 ]; then
    echo "Need at least 4 parameters! <BAM> <FRAGMENT> <CHROM_SIZES> <RESULT.bw>"
    exit 1
fi

BAM=$1
FRAGMENT=$2
CHROM_SIZES=$3
RESULT=$4

# Check configuration
[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util.sh

export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p "${TMPDIR}"

NAME=${BAM%%.bam}; NAME=${NAME##*/} # file name without extension

BDG=${RESULT/.bw/.bdg}
bash ${WASHU_ROOT}/scripts/bam2tags.sh ${BAM} ${FRAGMENT} |\
    awk -v OFS='\t' 'BEGIN{C="";S=0;E=0;X=0}
    {if(C!=$1||S!=$2||E!=$3){if(X!=0){print(C,S,E,X)};C=$1;S=$2;E=$3;X=1}else{X=X+1}}
    END{if(X!=0){print(C,S,E,X)}}' > ${BDG}
bash ${WASHU_ROOT}/scripts/bdg2bw.sh ${BDG} ${CHROM_SIZES}

# Cleanup
rm ${BDG}
# TMP dir cleanup:
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir