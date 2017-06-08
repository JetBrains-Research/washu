#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# CHPC (qsub) mock replacement
which qsub &>/dev/null || {
    qsub()
    {
        # LOAD args to $CMD
        while read -r line; do CMD+=$line; CMD+=$'\n'; done;
        QSUB_FILE=$(mktemp /tmp/qsub.XXXXXXXXXXXX.sh)
        echo "# This file was generated as QSUB MOCK" > $QSUB_FILE
        # MOCK for module command
        echo 'module() { echo "module $@"; } ' >> $QSUB_FILE
        echo "$CMD" >> $QSUB_FILE
        LOG=$(echo "$CMD" | grep "#PBS -o" | sed 's/#PBS -o //g')
        # Redirect both stderr and stdout to stdout then tee and then to stderr
        bash $QSUB_FILE 2>&1 | tee "$LOG" 1>&2
    }
}

# Small procedure to wait until all the tasks are finished on the qsub cluster
# Example of usage: wait_complete $TASKS, where $TASKS is a task ids returned by qsub.
wait_complete()
{
    echo "Waiting for tasks..."
    for TASK in $@
    do :
        echo -n "TASK: $TASK"
        # The task id is actually the first numbers in the string
        TASK_ID=$(echo ${TASK} | sed -e "s/\([0-9]*\).*/\1/")
        if [ ! -z "$TASK_ID" ]; then
            while qstat ${TASK_ID} &> /dev/null; do
                echo -n "."
                sleep 100
            done;
        fi
        echo
    done
    echo "Done."
}

# Checks for errors in logs, stops the world
check_logs()
{
    ERRORS=`find . -name "*.log" | xargs grep -i -e "err"`
    if [ ! -z "$ERRORS" ]; then
        echo "ERRORS found"
        echo "$ERRORS"
        exit 1
    fi
}

# Convert genome to macs2 species encoding
macs_species()
{
    GENOME=$1
    # Convert Genome build to macs2 species
    [[ ${GENOME} =~ ^hg[0-9]+$ ]] && SPECIES="hs"
    [[ ${GENOME} =~ ^mm[0-9]+$ ]] && SPECIES="mm"
    [[ -z "${SPECIES}" ]] && echo "Unknown species for macs: ${GENOME}" && exit 1
    echo "${SPECIES}"
}