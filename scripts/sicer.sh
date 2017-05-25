#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

which SICER.sh &>/dev/null || {
    echo "SICER not found! Download rseg: <http://home.gwu.edu/~wpeng/Software.htm>"
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
    INPUT=$(python $(dirname $0)/util.py find_input ${WORK_DIR}/${FILE})
    echo "${FILE} input: ${INPUT}"
    if [ ! -f "${INPUT}" ]; then
        echo "SICER requires control"
        continue
    fi

    NAME=${FILE%%.bam} # file name without extension


    # Create tmpfile in advance, because of interpolation of qsub call
    INPUT_FOLDER=$(mktemp -d)
    OUT_FOLDER=$(mktemp -d)
    FILE_TMP_BED=$(mktemp)

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
bedtools bamtobed -i ${FILE} | sort -k1,1 -k3,3n -k2,2n -k6,6 > ${FILE}.bed


# Use tmp files to reduced async problems with same input parallel processing
echo "${FILE}: control file found: ${INPUT}"
if [ ! -f ${INPUT}.bed ]; then
    bedtools bamtobed -i ${INPUT} | sort -k1,1 -k3,3n -k2,2n -k6,6 > ${FILE_TMP_BED}
    # Check that we are the first in async calls, not 100% safe
    if [ ! -f ${INPUT}.bed ]; then
        mv ${FILE_TMP_BED} ${INPUT}.bed
    fi
fi

# Ensure that folders exist
mkdir -p ${INPUT_FOLDER}
mkdir -p ${OUT_FOLDER}

cp ${FILE}.bed ${INPUT_FOLDER}
cp ${INPUT}.bed ${INPUT_FOLDER}
cd ${INPUT_FOLDER}

# Usage: SICER.sh [InputDir] [bed file] [control file] [OutputDir] [Species]
#   [redundancy threshold] [window size (bp)] [fragment size] [effective genome fraction] [gap size (bp)] [FDR]
# Defaults:
#   redundancy threshold    = 1
#   window size (bp)        = 200
#   fragment size           = 150
#   gap size (bp)           = 600

SICER.sh ${INPUT_FOLDER} ${FILE}.bed ${INPUT}.bed ${OUT_FOLDER} ${GENOME} 1 200 150 ${EFFECTIVE_GENOME_FRACTION} 600 ${FDR}
cp ${OUT_FOLDER}/island.bed ${WORK_DIR}/${NAME}_sicer_${FDR}.bed
ENDINPUT
)
    echo "FILE: ${FILE}; JOB: ${QSUB_ID}"
    TASKS="$TASKS $QSUB_ID"
done
wait_complete ${TASKS}
check_logs

# Cleanup
for FILE in $(find . -name '*.bam' | sed 's#./##g' | grep -v 'input')
do :
    INPUT=$(python $(dirname $0)/util.py find_input ${WORK_DIR}/${FILE})
    if [ -f "${INPUT}" ]; then
        rm ${INPUT.bed}
    fi
    rm ${FILE}.bed
done

echo "Done. Batch sicer: ${WORK_DIR} ${GENOME} ${CHROM_SIZES} ${FDR}"