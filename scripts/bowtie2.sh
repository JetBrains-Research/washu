#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Load technical stuff
source ~/work/washu/scripts/util.sh

if [ $# -lt 3 ]; then
    echo "Need 3 parameters! <WORK_DIR> <GENOME> <INDEXES>"
    exit 1
fi
WORK_DIR=$1
GENOME=$2
INDEXES=$3

echo "Batch Bowtie2: ${WORK_DIR} ${GENOME} ${INDEXES}"
cd ${WORK_DIR}

# Create soft link to indexes in working directory
INDEX_FILES=$(find ${INDEXES} -name "*.bt2*")
for F in ${INDEX_FILES[@]}; do TAG=${F##*/}; ln -s $F $TAG; done


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
#PBS -N bowtie2_${GENOME}_${NAME}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_bowtie2_${GENOME}.log

# Loading modules
module load bowtie2
module load samtools

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}
if [ -f "${FILE_PAIRED}" ]; then
    bowtie2 --sensitive-local -t -p 8 -S ${ID}.sam -x ${GENOME} -1 ${FILE} -2 ${FILE_PAIRED} &> ${ID}.bowtie2.log
else
    bowtie2 --sensitive-local -t -p 8 -S ${ID}.sam -x ${GENOME} ${FILE} &> ${ID}.bowtie2.log
fi
samtools view -bS -q 10 ${ID}.sam -o ${ID}_not_sorted.bam
samtools sort -@ 4 ${ID}_not_sorted.bam -o ${ID}.with_dup.bam
samtools rmdup ${ID}.with_dup.bam ${ID}.bam


# Cleanup
rm ${ID}.sam ${ID}_not_sorted.bam ${ID}.with_dup.bam

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

# Cleanup indexes soft link
rm *.bt2*

echo "Done. Batch Bowtie2: ${WORK_DIR} ${GENOME} ${INDEXES}"
