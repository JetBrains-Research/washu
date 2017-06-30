#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# MOCK for module command
which module &>/dev/null ||
    module() { echo "module $@"; }

# CHPC (qsub) mock replacement
which qsub &>/dev/null || {
    qsub()
    {
        # LOAD args to $CMD
        while read -r line; do CMD+=$line; CMD+=$'\n'; done;
        # MacOS cannot handle XXXX template with ".sh" suffix, also --suffix
        # option not available in BSD mktemp, so let's do some hack
        QSUB_FILE_PREFIX=$(mktemp "${TMPDIR:-/tmp/}qsub.XXXXXXXXXXXX")
        QSUB_FILE="${QSUB_FILE_PREFIX}.sh"
        rm ${QSUB_FILE_PREFIX}

        echo "# This file was generated as QSUB MOCK" > $QSUB_FILE
        # MOCK for module command
        echo 'module() { echo "module $@"; } ' >> $QSUB_FILE
        echo "$CMD" >> $QSUB_FILE
        LOG=$(echo "$CMD" | grep "#PBS -o" | sed 's/#PBS -o //g')
        # Redirect both stderr and stdout to stdout then tee and then to stderr
        bash $QSUB_FILE 2>&1 | tee "$LOG" 1>&2
    }
}

if which qsub &>/dev/null; then
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

    run_parallel()
    {
        # LOAD args to $CMD
        CMD=""
        while read -r line; do CMD+=$line; CMD+=$'\n'; done;

        QSUB_ID=$(qsub <<< "$CMD")
    }
else
    wait_complete()
    {
        echo "Waiting for tasks..."
        wait
        echo "Done."
    }

    run_parallel()
    {
        # Wait until less then 8 tasks running
        while [ `jobs | wc -l` -ge 8 ] ; do sleep 1 ; done

        CMD=""
        while read -r line; do CMD+=$line; CMD+=$'\n'; done;
        # MacOS cannot handle XXXX template with ".sh" suffix, also --suffix
        # option not available in BSD mktemp, so let's do some hack
        QSUB_FILE_PREFIX=$(mktemp "${TMPDIR:-/tmp/}qsub.XXXXXXXXXXXX")
        QSUB_FILE="${QSUB_FILE_PREFIX}.sh"
        rm ${QSUB_FILE_PREFIX}

        echo "# This file was generated as QSUB MOCK" > $QSUB_FILE
        # MOCK for module command
        echo 'module() { echo "module $@"; } ' >> $QSUB_FILE
        echo "$CMD" >> $QSUB_FILE
        LOG=$(echo "$CMD" | grep "#PBS -o" | sed 's/#PBS -o //g')

        # Redirect stdout and error to log file
        bash $QSUB_FILE &>"$LOG" &
        QSUB_ID="Task: ${QSUB_FILE} PID: $!"
    }
fi

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