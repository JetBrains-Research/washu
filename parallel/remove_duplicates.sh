#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Check configuration
[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util.sh

>&2 echo "Batch remove_duplicates $@"
if [[ $# -lt 1 ]]; then
    echo "Need at least 2 parameters! <PICARD_TOOLS_JAR> <WORK_DIR> [<WORK_DIR>]*"
    exit 1
fi

if [[ -n ${IS_TEST} ]]; then
    export JVM_MEM="-Xmx1g"
else
    export JVM_MEM="-Xmx12g"
fi

# Check Picard tools
PICARD_TOOLS_JAR=$1
if [[ ! -f "${PICARD_TOOLS_JAR}" ]]; then
    echo "Picard tools not found! Download Picard: <http://broadinstitute.github.io/picard/>"; exit 1;
fi
WORK_DIRS=${@:2}

PROCESSED=""
TASKS=()

for WORK_DIR in ${WORK_DIRS}; do
    WORK_DIR_NAME=${WORK_DIR##*/}
    cd ${WORK_DIR}

    for FILE in $(find . -name '*.bam' | grep -v _unique | sed 's#\./##g')
    do :
        NAME=${FILE%%.bam}
        UNIQUE_BAM=${NAME}_unique.bam
        METRICS=${NAME}_metrics.txt
        if [[ ! -f ${UNIQUE_BAM} ]]; then
            # Submit task
            run_parallel << SCRIPT
#!/bin/sh
#PBS -N unique_${WORK_DIR_NAME}_${NAME}
#PBS -l nodes=1:ppn=4,walltime=24:00:00,vmem=32gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_unique.log

cd ${WORK_DIR}
module load java

# PROBLEM: vmem is much bigger, however we face with the problem with bigger values:
# There is insufficient memory for the Java Runtime Environment to continue.
export _JAVA_OPTIONS="${JVM_MEM}"
source ${WASHU_ROOT}/parallel/util.sh
export TMPDIR=\$(type job_tmp_dir &>/dev/null && echo "\$(job_tmp_dir)" || echo "/tmp")

java -Djava.io.tmpdir="\${TMPDIR}" -jar ${PICARD_TOOLS_JAR} MarkDuplicates REMOVE_DUPLICATES=true INPUT=${FILE} OUTPUT=${UNIQUE_BAM} M=${METRICS}
SCRIPT
            echo "FILE: ${WORK_DIR_NAME}/${FILE}; TASK: ${QSUB_ID}"
            TASKS+=("$QSUB_ID")
        fi
    done
done
wait_complete ${TASKS[@]}
check_logs

>&2 echo "Done. Batch remove_duplicates $@"