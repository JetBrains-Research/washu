#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Check configuration
[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util.sh

>&2 echo "Batch SPAN $@"
if [ $# -lt 5 ]; then
    echo "Need >= 5 parameters! <SPAN_JAR_PATH> <WORK_DIR> <GENOME> <CHROM_SIZES> <Q> [<OUTPUT_DIR> [<GAP>]]"
    exit 1
fi

SPAN_JAR_PATH=$1
if [[ ! -f "${SPAN_JAR_PATH}" ]]; then
    >&2 echo "SPAN not found! Download SPAN: <http://artyomovlab.wustl.edu/aging/span.html>"; exit 1;
fi
WORK_DIR=$2
GENOME=$3
CHROM_SIZES=$4
Q=$5
OUTPUT_DIR=$6
if [ ! -d "$OUTPUT_DIR" ]; then
    OUTPUT_DIR=${WORK_DIR}
fi
GAP=$7
if [ -z "$GAP" ]; then
    GAP=5
fi

if [ -z "JAVA_OPTIONS" ]; then
    JAVA_OPTIONS="-Xmx8G"
fi

cd ${WORK_DIR}

TASKS=()
for FILE in $(find . -name '*.bam' | sed 's#\./##g' | grep -v 'input')
do :
    INPUT=$(python ${WASHU_ROOT}/scripts/util.py find_input ${WORK_DIR}/${FILE})
    echo "${FILE}: control file: ${INPUT}"

    NAME=${FILE%%.bam} # file name without extension
    ID=${NAME}_${Q}_${GAP}

    if [ ! -f ${ID}_peaks.bed ]; then
        # Submit task
        run_parallel << SCRIPT
#!/bin/sh
#PBS -N span_${ID}
#PBS -l nodes=1:ppn=4,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_span_${GENOME}.log

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

module load java
# PROBLEM: vmem is much bigger, however we face with the problem with bigger values:
# There is insufficient memory for the Java Runtime Environment to continue.
export _JAVA_OPTIONS="${JAVA_OPTIONS}"

if [ -f "${INPUT}" ]; then
    echo "${FILE}: control file found: ${INPUT}"
    java -jar ${SPAN_JAR_PATH} analyze -t ${FILE} -c ${INPUT} --chrom.sizes ${CHROM_SIZES} \
        --fdr ${Q} --gap ${GAP} \
        --peaks ${ID}_peaks.bed \
        --workdir ${OUTPUT_DIR} \
        --threads 4
else
    echo "${FILE}: no control file"
    java -jar ${SPAN_JAR_PATH} analyze -t ${FILE} --chrom.sizes ${CHROM_SIZES} \
        --fdr ${Q} --gap ${GAP} \
        --peaks ${ID}_peaks.bed \
        --workdir ${OUTPUT_DIR} \
        --threads 4
fi
SCRIPT

        echo "FILE: ${FILE}; TASK: ${QSUB_ID}"
        TASKS+=("$QSUB_ID")
    fi
done
wait_complete ${TASKS[@]}
check_logs

>&2 echo "Done. Batch SPAN $@"