#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Check Picard tools
if [[ ! -f ~/picard.jar ]]; then
    echo "Picard tools not found! Download Picard: <http://broadinstitute.github.io/picard/>"; exit 1;
fi

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

echo "Batch Bowtie2: ${WORK_DIR} ${GENOME} ${INDEXES} ${TRIM5}"
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
#PBS -l nodes=1:ppn=4,walltime=24:00:00,vmem=32gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_bowtie2_${GENOME}.log

# Loading modules
module load bowtie2
module load samtools

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

# Bowtie2 command line options
# -p/--threads <int> number of alignment threads to launch (1)
# -5/--trim5 <int>   trim <int> bases from 5'/left end of reads (0)
# --sensitive            -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 (default)

if [ -f "${FILE_PAIRED}" ]; then
    bowtie2 -p 4 --trim5 ${TRIM5} -S ${ID}.sam -x ${GENOME} -1 ${FILE} -2 ${FILE_PAIRED}
else
    bowtie2 -p 4 --trim5 ${TRIM5} -S ${ID}.sam -x ${GENOME} ${FILE}
fi
samtools view -bS -q 10 ${ID}.sam -o ${ID}_not_sorted.bam
samtools sort -@ 4 ${ID}_not_sorted.bam -o ${ID}.bam

module load java
# PROBLEM: vmem is much bigger, however we face with the problem with bigger values:
# There is insufficient memory for the Java Runtime Environment to continue.
export _JAVA_OPTIONS="-Xmx12g"
java -jar ~/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
    INPUT=${ID}.bam OUTPUT=${ID}_unique.bam M=${ID}_metrics.txt

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

# Cleanup indexes soft link
rm *.bt2*

echo "Done. Batch Bowtie2: ${WORK_DIR} ${GENOME} ${INDEXES} ${TRIM5}"
