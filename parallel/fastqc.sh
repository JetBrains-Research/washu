#!/usr/bin/env bash

###########################################################################
# Batch fastq-dump & multiqc:
#    Accepts list on one or many fastq containing directories.
#    In each <WORK_DIR> script runs fastqc for all its fastq/fastq.gz files
#    (single or paired ended) saves results to <WORK_DIR>/fastqc directory.
###########################################################################
# author oleg.shpynov@jetbrains.com
# author roman.chernyatchik@jetbrains.com

which multiqc &>/dev/null || {
    echo "fastq-dump not found! You can install it using:"
    echo "  conda install -c bioconda multiqc"
    echo "For further details see http://multiqc.info"
    exit 1
}

# Load technical stuff
source $(dirname $0)/../parallel/util.sh

>&2 echo "Batch fastqc $@"
if [ $# -lt 1 ]; then
    echo "Need at least one parameter! <WORK_DIR>"
    exit 1
fi
WORK_DIRS="$@"

for WORK_DIR in $(find . -type d -name "H*"); do :
    WORK_DIR_NAME=${WORK_DIR##*/}
    echo  "!${WORK_DIR_NAME}"

    cd ${WORK_DIR}
    find . -name "*_peaks.bed" -print0 | xargs -0 bash /home/user/work/washu/bed/union.sh | grep '|' | awk -v OFS='\t' '{print $1,$2,$3}' > ${WORK_DIR_NAME}_zinbra_weak_consensus.bed
    find . -name "*OD*_peaks.bed" -print0 | xargs -0 bash /home/user/work/washu/bed/union.sh | grep '|' | awk -v OFS='\t' '{print $1,$2,$3}' > ${WORK_DIR_NAME}_zinbra_ODS_weak_consensus.bed
    find . -name "*YD*_peaks.bed" -print0 | xargs -0 bash /home/user/work/washu/bed/union.sh | grep '|' | awk -v OFS='\t' '{print $1,$2,$3}' > ${WORK_DIR_NAME}_zinbra_YDS_weak_consensus.bed
done


TASKS=""
for WORK_DIR in ${WORK_DIRS}; do :
    cd ${WORK_DIR}
    WORK_DIR_NAME=${WORK_DIR##*/}

    mkdir -p "${WORK_DIR}/fastqc"

    for FILE in $(find . -regextype posix-extended -regex '.*\.f.*q(\.gz)?')
    do :
        FILE_NAME=${FILE##*/}
        NAME=${FILE_NAME%%.fast*} # file name without extension

        RESULTS_DIR="${WORK_DIR}/fastqc"
        if [ -d "${RESULTS_DIR}" ]; then
            echo "   [Skipped] ${RESULTS_DIR} was already processed"
            continue
        fi

        # Submit task
        run_parallel << SCRIPT
#!/bin/sh
#PBS -N fastqc_${WORK_DIR_NAME}_${NAME}
#PBS -l nodes=1:ppn=1,walltime=2:00:00,vmem=4gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_fastqc.log

# Loading modules
module load fastqc

# This is necessary because qsub default working dir is user home
# and our FILE is relative
cd ${WORK_DIR}

# Options:
# -o --outdir     Create all output files in the specified output directory.
#                     Please note that this directory must exist as the program
#                     will not create it. If this option is not set then the
#                      output file for each sequence file is created in the same
#                     directory as the sequence file which was processed.

# TODO: maybe use a couple of threads instead one?
# -t --threads    Specifies the number of files which can be processed
#                     simultaneously. Each thread will be allocated 250MB of
#                     memory so you shouldn't run more threads than your
#                     available memory will cope with, and not more than
#                      6 threads on a 32 bit machine
#
fastqc --outdir "${RESULTS_DIR}" "${FILE}"
SCRIPT
        echo "FILE: ${WORK_DIR_NAME}/${FILE}; TASK: ${QSUB_ID}"
        TASKS="$TASKS $QSUB_ID"
    done
done

wait_complete ${TASKS}
check_logs

for WORK_DIR in ${WORK_DIRS}; do :
    cd ${WORK_DIR}

    echo "Processing multiqc for: ${WORK_DIR}"
    #Options:
    # -f, --force           Overwrite any existing reports
    # -s, --fullnames       Do not clean the sample names (leave as full file name)
    # -o, --outdir TEXT     Create report in the specified output directory.
    multiqc -f -o "${WORK_DIR}" "${WORK_DIR}/fastqc"
done
>&2 echo "Done. Batch fastqc $@"
