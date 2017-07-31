#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

PICARD_TOOLS_JAR=~/picard.jar
# Check Picard tools
if [[ ! -f "${PICARD_TOOLS_JAR}" ]]; then
    echo "Picard tools not found! Download Picard: <http://broadinstitute.github.io/picard/>"; exit 1;
fi

# Load technical stuff, not available in qsub emulation
if [ -f "$(dirname $0)/util.sh" ]; then
    source "$(dirname $0)/util.sh"
fi

if [ $# -lt 1 ]; then
    echo "Need 1 parameters! <WORK_DIR>"
    exit 1
fi

WORK_DIR=$1

echo "Batch remove duplicates: ${WORK_DIR}"
cd ${WORK_DIR}

PROCESSED=""
TASKS=""

for FILE in $(find . -name '*.bam' | sed 's#\./##g')
do :
    NAME=${FILE%%.bam}
    UNIQUE_BAM=${NAME}_unique.bam
    METRICS=${NAME}_metrics.txt

    # Submit task
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N unique_${NAME}
#PBS -l nodes=1:ppn=4,walltime=24:00:00,vmem=32gb
#PBS -j oe
#PBS -o ${WORK_DIR}/${NAME}_unique.log

cd ${WORK_DIR}
module load java

# PROBLEM: vmem is much bigger, however we face with the problem with bigger values:
# There is insufficient memory for the Java Runtime Environment to continue.
export _JAVA_OPTIONS="-Xmx12g"
java -jar ${PICARD_TOOLS_JAR} MarkDuplicates REMOVE_DUPLICATES=true INPUT=${FILE} OUTPUT=${UNIQUE_BAM} M=${METRICS}

ENDINPUT
)
    echo "FILE: ${FILE}; TASK: ${QSUB_ID}"
    TASKS="$TASKS $QSUB_ID"
done
wait_complete ${TASKS}
check_logs

echo "Done. Batch remove duplicates: ${WORK_DIR}"