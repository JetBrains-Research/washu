#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

which SICER.sh &>/dev/null || {
    echo "SICER not found! Download SICER: <http://home.gwu.edu/~wpeng/Software.htm>"
    echo "Please refer to README for installation instructions, modify scripts, i.e."
    echo "sed -i 's#/home/data/SICER1.1#<YOUR_INSTALLATION_FOLDER>#g' SICER.sh"
    echo "SICER is python2 library, force it!"
    echo "sed -i 's#python#python2#g' SICER.sh"
    exit 1
}

# Load technical stuff, not available in qsub emulation
if [ -f "$(dirname $0)/util.sh" ]; then
    source "$(dirname $0)/util.sh"
fi

if [ $# -lt 4 ]; then
    echo "Need 4 parameters! <work_dir> <genome> <chrom.sizes> <FDR>"
    exit 1
fi

WORK_DIR=$1
GENOME=$2
CHROM_SIZES=$3
FDR=$4

echo "Batch sicer: ${WORK_DIR} ${GENOME} ${CHROM_SIZES} ${FDR}"
cd ${WORK_DIR}

EFFECTIVE_GENOME_FRACTION=$(python $(dirname $0)/util.py effective_genome_fraction ${GENOME} ${CHROM_SIZES})
echo "EFFECTIVE_GENOME_FRACTION: ${EFFECTIVE_GENOME_FRACTION}"

TASKS=""
for FILE in $(find . -name '*.bam' | sed 's#./##g' | grep -v 'input')
do :
    NAME=${FILE%%.bam} # file name without extension
    FILE_BED=${NAME}.bed

    INPUT=$(python $(dirname $0)/util.py find_input ${WORK_DIR}/${FILE})
    echo "${FILE} input: ${INPUT}"
    if [ ! -f "${INPUT}" ]; then
        echo "SICER requires control"
        continue
    fi
    INPUT_BED=${INPUT%%.bam}.bed


    # Create tmpfile in advance, because of interpolation of qsub call
    INPUT_FOLDER=${WORK_DIR}/tmp/${NAME}
    OUT_FOLDER=${WORK_DIR}/${NAME}/out
    # Create folders
    mkdir -p ${INPUT_FOLDER}
    mkdir -p ${OUT_FOLDER}

    # Submit task
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N sicer_${NAME}_${FDR}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_${FDR}_sicer.log

module load bedtools2

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

# SICER works with BED only
export LC_ALL=C
bedtools bamtobed -i ${FILE} | sort -k1,1 -k3,3n -k2,2n -k6,6 > ${INPUT_FOLDER}/${FILE_BED}

# Use tmp files to reduced async problems with same input parallel processing
echo "${FILE}: control file found: ${INPUT}"
if [ ! -f ${INPUT_BED} ]; then
    bedtools bamtobed -i ${INPUT} | sort -k1,1 -k3,3n -k2,2n -k6,6 > ${INPUT_FOLDER}/${INPUT_BED}
    # Check that we are the first in async calls, not 100% safe
    if [ ! -f ${INPUT_BED} ]; then
        cp ${INPUT_FOLDER}/${INPUT_BED} ${WORK_DIR}
    fi
else
    cp ${INPUT_BED} ${INPUT_FOLDER}
fi

cd ${INPUT_FOLDER}

# Usage: SICER.sh [InputDir] [bed file] [control file] [OutputDir] [Species]
#   [redundancy threshold] [window size (bp)] [fragment size] [effective genome fraction] [gap size (bp)] [FDR]
# Defaults:
#   redundancy threshold    = 1
#   window size (bp)        = 200
#   fragment size           = 150
#   gap size (bp)           = 600

SICER.sh ${INPUT_FOLDER} ${FILE_BED} ${INPUT_BED} ${OUT_FOLDER} ${GENOME} 1 200 150 ${EFFECTIVE_GENOME_FRACTION} 600 ${FDR}
cp -f ${OUT_FOLDER}/* ${WORK_DIR}
ENDINPUT
)
    echo "FILE: ${FILE}; JOB: ${QSUB_ID}"
    TASKS="$TASKS $QSUB_ID"
done
wait_complete ${TASKS}
check_logs

# Cleanup
rm -r ${WORK_DIR}/tmp
echo "Done. Batch sicer: ${WORK_DIR} ${GENOME} ${CHROM_SIZES} ${FDR}"