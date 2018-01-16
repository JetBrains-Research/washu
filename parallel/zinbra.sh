#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Load technical stuff
source $(dirname $0)/../parallel/util/util.sh
SCRIPT_DIR="$(project_root_dir)"

>&2 echo "Batch zinbra $@"
if [ $# -lt 5 ]; then
    echo "Need >= 5 parameters! <ZINBRA_JAR_PATH> <WORK_DIR> <GENOME> <CHROM_SIZES> <Q> [<OUTPUT_DIR> [<GAP>]]"
    exit 1
fi

ZINBRA_JAR_PATH=$1
if [[ ! -f "${ZINBRA_JAR_PATH}" ]]; then
    >&2 echo "ZINBRA not found! Download ZINBRA: <https://github.com/JetBrains-Research/zinbra>"; exit 1;
fi
WORK_DIR=$2
GENOME=$3
CHROM_SIZES=$4
Q=$5
OUTPUT_DIR=$6
if [ ! -d "$OUTPUT_DIR" ]; then
    OUTPUT_DIR=$WORK_DIR
fi

if [ $# -lt 7 ]; then
    GAP=""
else
    GAP="--gap $7"
fi

cd ${WORK_DIR}

TASKS=()
for FILE in $(find . -name '*.bam' | sed 's#\./##g' | grep -v 'input')
do :
    INPUT=$(python ${SCRIPT_DIR}/scripts/util.py find_input ${WORK_DIR}/${FILE})
    echo "${FILE}: control file: ${INPUT}"

    NAME=${FILE%%.bam} # file name without extension
    ID=${NAME}_${GENOME}_${Q}

    if [ ! -f ${ID}_peaks.bed ]; then
        # Submit task
        run_parallel << SCRIPT
#!/bin/sh
#PBS -N zinbra_${ID}
#PBS -l nodes=1:ppn=4,walltime=24:00:00,vmem=64gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_zinbra_${GENOME}.log

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

module load java
# PROBLEM: vmem is much bigger, however we face with the problem with bigger values:
# There is insufficient memory for the Java Runtime Environment to continue.
export _JAVA_OPTIONS="-Xmx30g"

if [ -f "${INPUT}" ]; then
    echo "${FILE}: control file found: ${INPUT}"
    java -cp ${ZINBRA_JAR_PATH} org.jetbrains.bio.zinbra.ZinbraCLA analyze -t ${FILE} -c ${INPUT} \
        --chrom.sizes ${CHROM_SIZES} --fdr ${Q} --bed ${ID}_peaks.bed \
        --output ${OUTPUT_DIR} \
        --threads=4 ${GAP}

else
    echo "${FILE}: no control file"
    java -cp ${ZINBRA_JAR_PATH} org.jetbrains.bio.zinbra.ZinbraCLA analyze -t ${FILE} \
        --chrom.sizes ${CHROM_SIZES} --fdr ${Q} --bed ${ID}_peaks.bed \
        --output ${OUTPUT_DIR} \
        --threads=4 ${GAP}
fi

module load bedtools2
# Compute Reads in Peaks
bash ${SCRIPT_DIR}/scripts/rip.sh ${FILE} ${ID}_peaks.bed
SCRIPT

        echo "FILE: ${FILE}; TASK: ${QSUB_ID}"
        TASKS+=("$QSUB_ID")
    fi
done
wait_complete ${TASKS[@]}
check_logs

>&2 echo "Done. Batch zinbra $@"
