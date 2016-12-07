#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Load technical stuff
source ~/work/washu/scripts/util.sh

if [ $# -lt 4 ]; then
    echo "Need 4 parameters! <WORK_DIR> <GENOME> <INDEXES> <TRIM5>"
    exit 1
fi
WORK_DIR=$1
GENOME=$2
INDEXES=$3
TRIM5=$4

echo "Batch Bowtie: ${WORK_DIR} ${GENOME} ${INDEXES} ${TRIM5}"
cd ${WORK_DIR}

PROCESSED=""
TASKS=""
for FILE in $(find . -name '*.f*q' -printf '%P\n' | sort)
do :
    if $(echo "${PROCESSED[@]}"  | fgrep -q "${FILE}");
    then
        echo "$FILE was already processed"
        continue
    fi

    # Assumption: the only difference between paired-end read files is _1 and _2
    FILE_PAIRED=""
    if $(echo "${FILE_PAIRED[@]}"  | fgrep -q "_1");
    then
        PREFIX=${FILE%%_1.*}
        SUFFIX=${FILE##*_1}
        FILE_PAIRED="${PREFIX}_2${SUFFIX}"
    fi

    # Setup correct name
    if [ -f "${FILE_PAIRED}" ]; then
        echo "PAIRED END reads detected: $FILE and $FILE_PAIRED"
        # Mark it as already processed
        PROCESSED="${PROCESSED} ${FILE_PAIRED}"
        NAME="${PREFIX}${SUFFIX}"
    else
        NAME=${FILE%%.f*q}
    fi
    ID=${NAME}_${GENOME}

    # Submit task
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N bowtie_${GENOME}_${NAME}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_bowtie_${GENOME}.log

# Loading modules
module load bowtie
module load samtools

export BOWTIE_INDEXES=${INDEXES}
# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}
if [ -f "${FILE_PAIRED}" ]; then
    bowtie -p 8 -St -m 1 -v 3 --trim5 ${TRIM5} --best --strata ${GENOME} -1 ${FILE} -2 ${FILE_PAIRED} ${ID}.sam
else
    bowtie -p 8 -St -m 1 -v 3 --trim5 ${TRIM5} --best --strata ${GENOME} ${FILE} ${ID}.sam
fi
samtools view -bS ${ID}.sam -o ${ID}_not_sorted.bam
samtools sort ${ID}_not_sorted.bam -o ${ID}.bam

# Cleanup
rm ${ID}.sam ${ID}_not_sorted.bam

ENDINPUT
)
    if [ -f "${FILE_PAIRED}" ]; then
        echo "FILE: ${FILE} PAIRED ${FILE_PAIRED}; JOB: ${QSUB_ID}"
    else
        echo "FILE: ${FILE}; JOB: ${QSUB_ID}"
    fi
    TASKS="$TASKS $QSUB_ID"
done
wait_complete ${TASKS}
check_logs

echo "Done. Batch Bowtie: ${WORK_DIR} ${GENOME} ${INDEXES} ${TRIM5}"