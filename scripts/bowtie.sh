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

echo "Batch Bowtie: ${WORK_DIR} ${GENOME} ${INDEXES} ${TRIM5}"
cd ${WORK_DIR}

PROCESSED=""
TASKS=""

# Fails with large indexes, create soft link to indexes in working directory as a workaround
# export BOWTIE_INDEXES=${INDEXES}
if [ ! -z "${WORK_DIR}/indexes" ]; then
    ln -s ${INDEXES} ${WORK_DIR}/indexes
fi

# Bowtie fails with large indexes, explicitly set
if [ -f indexes/${GENOME}.1.ebwtl ]; then
    INDEX_ARG="--large-index"
else
    INDEX_ARG=""
fi

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
#PBS -l nodes=1:ppn=4,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_bowtie_${GENOME}.log

# Loading modules
module load bowtie
module load samtools

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

# Bowtie command line options used
# -p/--threads <int> number of alignment threads to launch (default: 1)
# -S/--sam           write hits in SAM format
# -t/--time          print wall-clock time taken by search phases
# -m <int>           suppress all alignments if > <int> exist (def: no limit)
# -v <int>           report end-to-end hits w/ <=v mismatches; ignore qualities
# -5/--trim5 <int>   trim <int> bases from 5' (left) end of reads
# --best             hits guaranteed best stratum; ties broken by quality
# --strata           hits in sub-optimal strata aren't reported (requires --best)

if [ -f "${FILE_PAIRED}" ]; then
    bowtie -p 4 -St -m 1 -v 3 --trim5 ${TRIM5} --best --strata ${INDEX_ARG} ${GENOME} -1 ${FILE} -2 ${FILE_PAIRED} ${ID}.sam
else
    bowtie -p 4 -St -m 1 -v 3 --trim5 ${TRIM5} --best --strata ${INDEX_ARG} ${GENOME} ${FILE} ${ID}.sam
fi
samtools view -bS ${ID}.sam -o ${ID}_not_sorted.bam
samtools sort ${ID}_not_sorted.bam -o ${ID}.bam

# Remove duplicated reads
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
if [ -z "${WORK_DIR}/indexes" ]; then
    rm ${WORK_DIR}/indexes
fi

echo "Done. Batch Bowtie: ${WORK_DIR} ${GENOME} ${INDEXES} ${TRIM5}"