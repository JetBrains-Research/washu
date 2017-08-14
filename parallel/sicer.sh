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

>&2 echo "Batch sicer $@"
if [ $# -lt 4 ]; then
    echo "Need at least 4 parameters! <work_dir> <genome> <chrom.sizes> <FDR> [window size (bp)] [fragment size] [gap size (bp)]"
    exit 1
fi

WINDOW_SIZE=200
FRAGMENT_SIZE=150
GAP_SIZE=600

if [ $# -ge 5 ]; then WINDOW_SIZE=$5 ; fi
if [ $# -ge 6 ]; then FRAGMENT_SIZE=$6 ; fi
if [ $# -ge 7 ]; then GAP_SIZE=$7 ; fi

WORK_DIR=$1
GENOME=$2
CHROM_SIZES=$3
FDR=$4

cd ${WORK_DIR}

EFFECTIVE_GENOME_FRACTION=$(python $(dirname $0)/../scripts/util.py effective_genome_fraction ${GENOME} ${CHROM_SIZES})
echo "EFFECTIVE_GENOME_FRACTION: ${EFFECTIVE_GENOME_FRACTION}"

if [ -z "${EFFECTIVE_GENOME_FRACTION}" ]; then
    echo "EFFECTIVE_GENOME_FRACTION is not determined"
    exit 1
fi

TASKS=""
for FILE in $(find . -name '*.bam' | sed 's#\./##g' | grep -v 'input')
do :
    NAME=${FILE%%.bam} # file name without extension
    FILE_BED=${NAME}.bed

    INPUT=$(python $(dirname $0)/../scripts/util.py find_input ${WORK_DIR}/${FILE})
    echo "${FILE} input: ${INPUT}"
    if [ ! -f "${INPUT}" ]; then
        echo "SICER requires control"
        continue
    fi
    INPUT_BED=${INPUT%%.bam}.bed


    # Create tmpfile in advance, because of interpolation of qsub call
    INPUT_FOLDER=${WORK_DIR}/sicer_tmp/${NAME}
    OUT_FOLDER=${INPUT_FOLDER}/out
    # Create folders
    mkdir -p ${INPUT_FOLDER}
    mkdir -p ${OUT_FOLDER}

    # Submit task
    run_parallel << ENDINPUT
#!/bin/sh
#PBS -N sicer_${NAME}_${FDR}
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=16gb
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


echo "SICER.sh ${INPUT_FOLDER} ${FILE_BED} ${INPUT_BED} ${OUT_FOLDER} ${GENOME} 1 ${WINDOW_SIZE} ${FRAGMENT_SIZE} ${EFFECTIVE_GENOME_FRACTION} ${GAP_SIZE} ${FDR}"

SICER.sh ${INPUT_FOLDER} ${FILE_BED} ${INPUT_BED} ${OUT_FOLDER} ${GENOME} 1 ${WINDOW_SIZE} ${FRAGMENT_SIZE} ${EFFECTIVE_GENOME_FRACTION} ${GAP_SIZE} ${FDR}
cp -f ${OUT_FOLDER}/* ${WORK_DIR}

# Output files:
# -normalized.wig. This file is in WIG format and could be uploaded to the UCSC genome browser for visualization of the ChIP library with redundancy-removed but before island-filtering (ES_H3K27me3-1-removed.bed) with desired window size.
# .scoreisland. This file stores all identified islands with respective scores and could be used to evaluate the choice of gap size.
# -islands-summary. This 8-column file is a summary of all islands identified in the ChIP library. The format is chromosome, start, end, read count in ChIP library, read count in control library, p-value, fold change, FDR).
# -island.bed. This file has 4 columns and contains information about all identified ChIP-enriched region information under filtering parameters (window size 200, gap size 600, FDR cutoff 1e-3). The 4 columns are chromosome, start, end, and read count in ChIP library.
# -islandfiltered.bed. This file is in BED format and contains all reads that are within significant islands.
# -islandfiltered-normalized.wig. This file is in WIG format and could be uploaded directly to UCSC genome browser for visualization of the island-filtered ChIP library.
# -removed.bed. This file contains reads after redundancy removed.

cd ${WORK_DIR}

# Compute Reads in Peaks
bash $(dirname $0)/../reports/rip.sh ${FILE} ${NAME}*island.bed
ENDINPUT

    echo "FILE: ${FILE}; TASK: ${QSUB_ID}"
    TASKS="$TASKS $QSUB_ID"
done
wait_complete ${TASKS}
check_logs

# Cleanup
rm -r ${WORK_DIR}/sicer_tmp
>&2 echo "Done. Batch sicer $@"