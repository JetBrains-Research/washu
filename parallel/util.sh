#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# MOCK for module command
type module &>/dev/null ||
    module() { echo "[mock] module $@"; }

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
        echo 'which module &>/dev/null || module() { echo "module $@"; }' >> $QSUB_FILE

        echo "$CMD" >> $QSUB_FILE
        LOG=$(echo "$CMD" | grep "#PBS -o" | sed 's/#PBS -o //g')

        # Log tasks script path to STDERR, STDOUT will be treated as qsub
        # return value which is supposed to be task name and is handled
        # by wait_complete
        echo "Running TASK: ${QSUB_FILE}" 1>&2
        # Redirect both stderr and stdout to stdout then tee and then to stderr
        bash $QSUB_FILE &> "$LOG"
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

    # Use function to get rid of command substitution.
    # Command substitution not works well with parallel execution.
    run_parallel()
    {
        # LOAD args to $CMD
        CMD=""
        while read -r line; do CMD+=$line; CMD+=$'\n'; done;

        # Return through global variable here, because we can't use command substitution.
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
        while [ $(jobs | wc -l) -ge 8 ] ; do sleep 1 ; done

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

# Convert path to absolute path and expand all symlinks
function expand_path() {
    # expand ".." and "." including trailing case
    # based on https://stackoverflow.com/questions/3915040/bash-fish-command-to-print-absolute-path-to-a-file
    TARGET_FILE="$(pwd)/$1"

    if [ -d "$1" ]; then
        # dir
        TARGET_FILE="$(cd "$1"; pwd)"
    elif [ -f "$1" ]; then
        # file
        if [[ $1 == */* ]]; then
            TARGET_FILE="$(cd "${1%/*}"; pwd)/${1##*/}"
        fi
    fi

    # resolve symlinks:
    # replacement for `readlink -f path` which isn't available in MacOS
    # http://stackoverflow.com/questions/1055671/how-can-i-get-the-behavior-of-gnus-readlink-f-on-a-mac
    cd "$(dirname ${TARGET_FILE})"
    TARGET_FILE=`basename ${TARGET_FILE}`

    # Iterate down a (possible) chain of symlinks
    while [ -L "$TARGET_FILE" ]
    do
        TARGET_FILE="$(readlink ${TARGET_FILE})"
        cd "$(dirname ${TARGET_FILE})"
        TARGET_FILE="$(basename ${TARGET_FILE})"
    done

    # Compute the canonicalized name by finding the physical path
    # for the directory we're in and appending the target file.
    PHYS_DIR="$(pwd -P)"
    echo "${PHYS_DIR}/${TARGET_FILE}"
}

PROJECT_ROOT_DIR="$(expand_path "$(dirname $0)/..")"
project_root_dir() {
    echo ${PROJECT_ROOT_DIR}
}

# Checks for errors in logs, stops the world
check_logs()
{
    ERRORS=`find . -name "*.log" | xargs grep -i -e "error"`
    if [ ! -z "$ERRORS" ]; then
        echo "ERRORS found"
        echo "$ERRORS"
    fi
}