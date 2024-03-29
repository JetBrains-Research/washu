#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Check configuration
[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util.sh

>&2 echo "Batch bowtie $@"
if [[ $# -lt 4 ]]; then
    echo "Need 4 parameters! <GENOME> <INDEXES> <TRIM5> <WORK_DIR> [<WORK_DIR>]"
    exit 1
fi

GENOME=$1
INDEXES=$2
TRIM5=$3
WORK_DIRS=${@:4}

# Configure TRIM_ARGS to avoid bowtie segfault in 0 case
if [[ ${TRIM5} -gt 0 ]]; then
    TRIM_ARGS="--trim5 ${TRIM5}"
else
    TRIM_ARGS=""
fi

TASKS=()
PROCESSED=()

for WORK_DIR in ${WORK_DIRS}; do :
    WORK_DIR_NAME=${WORK_DIR##*/}
    cd ${WORK_DIR}

    # Fails with large indexes, create soft link to indexes in working directory as a workaround
    # export BOWTIE_INDEXES=${INDEXES}
    if [[ ! -d "${WORK_DIR}/indexes" ]]; then
        ln -s ${INDEXES} ${WORK_DIR}/indexes
    fi

    # Bowtie fails with large indexes, explicitly set
    if [[ -f indexes/${GENOME}.1.ebwtl ]]; then
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


        # Assumption: the only difference between paired-end read files is _1 and _2 / _R1 and _R2
        FILE_PAIRED=""
        if $(echo "${FILE}"  | fgrep -q "_1"); then
            FILE_PAIRED=$(echo "${FILE}" | sed 's#_1#_2#g')
            NAME=$(echo $(basename ${FILE%%.f*q}) | sed 's#_1##g')
        elif $(echo "${FILE}" | fgrep -q "_R1"); then
            FILE_PAIRED=$(echo "${FILE}" | sed 's#_R1#_R2#g')
            NAME=$(echo $(basename ${FILE%%.f*q}) | sed 's#_R1##g')
        else
            NAME=$(basename ${FILE%%.f*q})
        fi

        # Check FILE_PAIRED
        if [[ -f "${FILE_PAIRED}" ]]; then
            echo "PAIRED END reads detected: $FILE and $FILE_PAIRED"
            # Mark it as already processed
            PROCESSED+=("${FILE_PAIRED}")
        fi

        ID=${NAME}_${GENOME}
        BAM_NAME="${ID}.bam"

        if [[ -f "${BAM_NAME}" ]]; then
            echo "   [Skipped] ${WORK_DIR}/${BAM_NAME} already exists."
            continue
        fi

        # Submit task
        run_parallel << SCRIPT
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

if [[ -f "${FILE_PAIRED}" ]]; then
    bowtie -p 4 -St -m 1 -v 3 ${TRIM_ARGS} --best --strata ${INDEX_ARG} ${GENOME} -1 ${FILE} -2 ${FILE_PAIRED} ${ID}.sam
else
    bowtie -p 4 -St -m 1 -v 3 ${TRIM_ARGS} --best --strata ${INDEX_ARG} ${GENOME} ${FILE} ${ID}.sam
fi
samtools view -bS ${ID}.sam -o ${ID}_not_sorted.bam
samtools sort ${ID}_not_sorted.bam -o ${BAM_NAME}

# Cleanup
rm ${ID}.sam ${ID}_not_sorted.bam

SCRIPT
        if [[ -f "${FILE_PAIRED}" ]]; then
            echo "FILE: ${WORK_DIR_NAME}/${FILE} PAIRED ${FILE_PAIRED}; TASK: ${QSUB_ID}"
        else
            echo "FILE: ${WORK_DIR_NAME}/${FILE}; TASK: ${QSUB_ID}"
        fi
        TASKS+=("$QSUB_ID")
    done
done
wait_complete ${TASKS[@]}
check_logs

# Cleanup indexes soft link
for WORK_DIR in ${WORK_DIRS}; do :
    if [[ -d "${WORK_DIR}/indexes" ]]; then
        rm ${WORK_DIR}/indexes
    fi
done

>&2 echo "Done. Batch bowtie $@"