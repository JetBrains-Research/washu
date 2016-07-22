#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

# Small procedure to wait until all the tasks are finished on the qsub cluster
# Example of usage: wait_complete $TASKS, where $TASKS is a task ids returned by qsub.
wait_complete()
{
    echo "Waiting for tasks."
    for JOB in $@
    do :
        echo -n "JOB: $JOB"
        # The job id is actually the first numbers in the string
        JOB_ID=$(echo ${JOB} | sed -e "s/\([0-9]*\).*/\1/")
        if [ ! -z "$JOB_ID" ]; then
            while qstat ${JOB_ID} &> /dev/null; do
                echo -n "."
                sleep 5
            done;
        fi
        echo
    done
    echo "Done."
}

# Checks for errors in logs, stops the world
check_logs()
{

    ERRORS=`find . -name "*.log" | xargs grep -i -e "err|"`
    if [ ! -z "$ERRORS" ]; then
        echo "ERRORS found"
        echo "$ERRORS"
        exit 1
    fi
}