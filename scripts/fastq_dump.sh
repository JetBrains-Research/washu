#!/usr/bin/env bash

###########################################################################
# Batch fastq-dump:
#    Accepts list on one or many SRA containing directories. In each <WORK_DIR>
#    script runs fastq-dump for all its *.sra files and saves results
#    to <WORK_DIR>/fastq/ directory. If file already exists it will be skipped
#
#    Ensure tat fastq-dump is installed or install it using:
#       conda install -c bioconda sra-tools
#    For further details see https://www.ncbi.nlm.nih.gov/books/NBK158900
###########################################################################
# author roman.chernyatchik@jetbrains.com

# Load technical stuff, not available in qsub emulation
if [ -f "$(dirname $0)/util.sh" ]; then
    source "$(dirname $0)/util.sh"
fi

if [ $# -lt 1 ]; then
    echo "Need at least one parameter! <WORK_DIR>"
    exit 1
fi

WORK_DIRS="$@"
TASKS=""

echo "Batch fastq-dump: ${WORK_DIRS}"

for WORK_DIR in ${WORK_DIRS}; do :
    cd ${WORK_DIR}

    OUTDIR="${WORK_DIR}/fastq"
    mkdir -p ${OUTDIR}

    WORK_DIR_NAME=${WORK_DIR##*/}
    for FILE in $(find ${WORK_DIR} -name '*.sra' | sort); do :
        FILE_NAME=${FILE##*/}
        FILE_NAME_WO_EXT=${FILE_NAME%.sra}

        FASTQ_FILE_PREFIX="${OUTDIR}/${FILE_NAME_WO_EXT}"
        if [ -f "${FASTQ_FILE_PREFIX}.fastq.gz" ] ||
         ( [ -f "${FASTQ_FILE_PREFIX}_1.fastq.gz" ] &&
         [ -f "${FASTQ_FILE_PREFIX}_2.fastq.gz" ]); then
            echo "  $FILE_NAME was already processed"
            continue
        fi

        # Submit task
        QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N fastq-dump_${WORK_DIR_NAME}_${FILE_NAME_WO_EXT}
#PBS -l nodes=1:ppn=4,walltime=24:00:00,vmem=32gb
#PBS -j oe
#PBS -o ${OUTDIR}/${FILE_NAME_WO_EXT}_fastq_dump.log

# Loading modules
module load sratoolkit

# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}

NLINES=\$(fastq-dump --maxSpotId 1 --split-spot --stdout ${FILE} 2>/dev/null | wc -l)
if [ \$NLINES -eq 4 ]; then
    echo "${FILE}: SE reads"
    SPLIT_FILES_OPTION=""
elif [ \$NLINES -eq 8 ]; then
    echo "${FILE}: PE reads"
    SPLIT_FILES_OPTION=" --split-files"
else
    echo "${FILE}: cannot detect whether single or paired reads"
    exit 1
fi

# fastq-dump command line options used
# -L|--log-level <level>  Logging level as number or enum string One of (fatal|sys|int|err|warn|info)
# -O|--outdir <path>      Output directory, default is working directory '.' )
# -B|--dumpbase           Formats sequence using base space (default for other than SOLiD).
# --gzip                  Compress output using gzip
# -O|--outdir <path>      Output directory, default is working

fastq-dump --log-level err --dumpbase --gzip --outdir ${OUTDIR}\${SPLIT_FILES_OPTION} --helicos ${FILE}

ENDINPUT
)
        echo "FILE: ${FILE}; TASK: ${QSUB_ID}"
        TASKS="$TASKS $QSUB_ID"
    done
done

wait_complete ${TASKS}
check_logs

echo "Done. Batch fastq-dump: ${WORK_DIRS}"