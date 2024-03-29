#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Check configuration
[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util.sh

>&2 echo "Batch bowtie2 $@"
if [[ $# -lt 4 ]]; then
    echo "Need at least 4 parameters! <GENOME> <INDEXES> <TRIM5> <WORK_DIR> [<WORK_DIR>]"
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

    # Create soft link to indexes in working directory
    INDEX_FILES=$(find ${INDEXES} -name "*.bt2*")
    for F in ${INDEX_FILES[@]}; do TAG=${F##*/}; ln -s $F $TAG; done

    for FILE in $(find . -regextype posix-extended -regex '.*\.f.*q(\.gz)?' | sort)
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
#PBS -N bowtie2_${GENOME}_${WORK_DIR_NAME}_${NAME}
#PBS -l nodes=1:ppn=4,walltime=24:00:00,vmem=32gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_bowtie2_${GENOME}.log

# Loading modules
module load bowtie2
module load samtools

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

# Bowtie2 command line options
# bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i>} [-S <sam>]
#
#  <bt2-idx>  Index filename prefix (minus trailing .X.bt2).
#             NOTE: Bowtie 1 and Bowtie 2 indexes are not compatible.
#  <m1>       Files with #1 mates, paired with files in <m2>.
#             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
#  <m2>       Files with #2 mates, paired with files in <m1>.
#             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
#  <r>        Files with unpaired reads.
#             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
#  <i>        Files with interleaved paired-end FASTQ reads
#             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
#  <sam>      File for SAM output (default: stdout)
#
# -p/--threads <int> number of alignment threads to launch (1)
# -5/--trim5 <int>   trim <int> bases from 5'/left end of reads (0)
# --sensitive            -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 (default)

if [[ -f "${FILE_PAIRED}" ]]; then
    bowtie2 -p 4 ${TRIM_ARGS} -S ${ID}.sam -x ${GENOME} -1 ${FILE} -2 ${FILE_PAIRED}
else
    bowtie2 -p 4 ${TRIM_ARGS} -S ${ID}.sam -x ${GENOME} -U ${FILE}
fi
samtools view -bS -q 10 ${ID}.sam -o ${ID}_not_sorted.bam
samtools sort -@ 4 ${ID}_not_sorted.bam -o ${BAM_NAME}

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
    cd ${WORK_DIR}
    rm *.bt2*
done

>&2 echo "Done. Batch bowtie2 $@"
