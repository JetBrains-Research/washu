#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Load technical stuff, not available in qsub emulation
if [ -f "$(dirname $0)/util.sh" ]; then
    source "$(dirname $0)/util.sh"
fi

>&2 echo "Batch bowtie $@"
if [ $# -lt 4 ]; then
    echo "Need 4 parameters! <GENOME> <INDEXES> <TRIM5> <WORK_DIR> [<WORK_DIR>]"
    exit 1
fi

GENOME=$1
INDEXES=$2
TRIM5=$3
WORK_DIRS=${@:4}

TASKS=""
PROCESSED=""

for WORK_DIR in ${WORK_DIRS}; do :
    WORK_DIR_NAME=${WORK_DIR##*/}
    cd ${WORK_DIR}

    # Fails with large indexes, create soft link to indexes in working directory as a workaround
    # export BOWTIE_INDEXES=${INDEXES}
    if [ ! -d "${WORK_DIR}/indexes" ]; then
        ln -s ${INDEXES} ${WORK_DIR}/indexes
    fi

    # Bowtie fails with large indexes, explicitly set
    if [ -f indexes/${GENOME}.1.ebwtl ]; then
        INDEX_ARG="--large-index"
    else
        INDEX_ARG=""
    fi

    for FILE in $(find . -name '*.f*q' | sed 's#\./##g' | sort)
    do :
        if $(echo "${PROCESSED[@]}"  | fgrep -q "${FILE}");
        then
            echo "$FILE was already processed"
            continue
        fi

        # Assumption: the only difference between paired-end read files is _1 and _2
        FILE_PAIRED=""
        if $(echo "${FILE}"  | fgrep -q "_1");
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
        NAME="${NAME##*/}" # if relative path - trim folders
        ID=${NAME}_${GENOME}

        echo "ID=${ID}"
        echo "wd=${WORK_DIR}"
        echo "fp=${FILE_PAIRED}"
        echo "file=${FILE}"

        # Submit task
        QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N bowtie_${GENOME}_${WORK_DIR_NAME}_${NAME}
#PBS -l nodes=1:ppn=4,walltime=24:00:00,vmem=32gb
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

# Cleanup
rm ${ID}.sam ${ID}_not_sorted.bam

ENDINPUT
)
        if [ -f "${FILE_PAIRED}" ]; then
            echo "FILE: ${WORK_DIR_NAME}/${FILE} PAIRED ${FILE_PAIRED}; TASK: ${QSUB_ID}"
        else
            echo "FILE: ${WORK_DIR_NAME}/${FILE}; TASK: ${QSUB_ID}"
        fi
        TASKS="$TASKS $QSUB_ID"
    done
done
wait_complete ${TASKS}
check_logs

# Cleanup indexes soft link
for WORK_DIR in ${WORK_DIRS}; do :
    if [ -d "${WORK_DIR}/indexes" ]; then
        rm ${WORK_DIR}/indexes
    fi
done

>&2 echo "Done. Batch bowtie $@"