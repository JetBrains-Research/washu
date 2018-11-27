#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

which rsem-prepare-reference &>/dev/null || { echo "RSEM not found! Download RSEM: <https://github.com/deweylab/RSEM>"; exit 1; }

# Check configuration
[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util.sh

>&2 echo "index-rsem $@"
if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <GENOME> <FOLDER>"
    exit 1
fi
GENOME=$1
FOLDER=$2

RSEM_INDEX_FOLDER="${FOLDER}/rsem"
echo "RSEM index folder: ${RSEM_INDEX_FOLDER}"
if [ -d "${RSEM_INDEX_FOLDER}" ]; then
    echo "Indexes already exist at ${RSEM_INDEX_FOLDER}"
    exit 0
fi

STAR_INDEX_FOLDER="${FOLDER}/star"
echo "STAR index folder: ${STAR_INDEX_FOLDER}"
if [ ! -d "${STAR_INDEX_FOLDER}" ]; then
    echo "STAR Indexes not found ${STAR_INDEX_FOLDER}"
    exit 1
fi

mkdir -p ${RSEM_INDEX_FOLDER}
run_parallel << SCRIPT
#!/bin/sh
#PBS -N rsem_indexes_${GENOME}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=64gb
#PBS -j oe
#PBS -o ${FOLDER}/${GENOME}_rsem_indexes.log

# This is necessary because qsub default working dir is user home
cd ${RSEM_INDEX_FOLDER}

rsem-prepare-reference --gtf ${STAR_INDEX_FOLDER}/${GENOME}.gtf ${FOLDER} ${GENOME}
SCRIPT

wait_complete ${QSUB_ID}
check_logs

>&2 echo "Done. index-rsem $@"