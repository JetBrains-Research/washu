#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Load technical stuff
source ~/work/washu/scripts/util.sh

WORK_DIR=$1
GENOME=$2
Q=$3

echo "Batch macs2 broad: ${WORK_DIR} ${GENOME} ${Q}"

# Convert Genome build to macs2 species
[[ ${GENOME} =~ ^hg[0-9]+$ ]] && SPECIES="hs"
[[ ${GENOME} =~ ^mm[0-9]+$ ]] && SPECIES="mm"
[[ -z "$SPECIES" ]] && echo "Unknown species for macs: $GENOME" && exit 1

cd ${WORK_DIR}

TASKS=""
for FILE in $(find . -type f -name '*.bam' -printf '%P\n')
do :
    # Convention over configuration: we assume that input has the same naming scheme as chromatin marks
    if [[ ! ${FILE} =~ ^.*input.*$ ]]; then
        INPUT=$(echo ${FILE} | sed -e 's/k27me3/input/g' -e 's/k27ac/input/g' -e 's/k4me3/input/g')
    else
        INPUT="" # No input for itself
    fi
    echo "${FILE} input: ${INPUT}"

    NAME=${FILE%%.bam} # file name without extension
    ID=${NAME}_broad_${Q}

    # Submit task
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N macs2_${ID}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${ID}_macs2_broad_${GENOME}.log

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}
if [ -f "${INPUT}" ]; then
    echo "Input file found: ${INPUT}"
    /home/oshpynov/miniconda2/bin/macs2 callpeak -t ${FILE} -c ${INPUT} -f BAM -g ${SPECIES} -n ${ID} -B --broad --broad-cutoff ${Q}

else
    echo "No input file"
    /home/oshpynov/miniconda2/bin/macs2 callpeak -t ${FILE} -f BAM -g ${SPECIES} -n ${ID} -B --broad --broad-cutoff ${Q}
fi
ENDINPUT
)
    echo "$FILE: $QSUB_ID"
    TASKS="$TASKS $QSUB_ID"
done
wait_complete ${TASKS}
check_logs

echo "Batch macs2 broad: ${WORK_DIR} ${GENOME} ${Q}"