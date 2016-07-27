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
    ERRORS=`find . -name "*.log" | xargs grep -i -e "err"`
    if [ ! -z "$ERRORS" ]; then
        echo "ERRORS found"
        echo "$ERRORS"
        exit 1
    fi
}

# Find control/input for given file
macs2_find_control()
{
    FILE=$1
    # Convention over configuration: we assume that input has the same naming scheme as chromatin marks
    if [[ ! ${FILE} =~ ^.*input.*$ ]]; then
        DONOR=$(echo ${FILE} | sed -e "s/.*\(donor[0-9]\).*/\1/")
        echo "${FILE} donor: ${DONOR}"
        INPUTS=$(find . -name "*${DONOR}_input*.bam" -printf '%P\n')
        INPUT=${INPUTS[0]}
    else
        INPUT="" # No input for itself
    fi
    echo "${INPUT}"
}

# Convert genome to macs2 species encoding
macs2_species()
{
    GENOME=$1
    # Convert Genome build to macs2 species
    [[ ${GENOME} =~ ^hg[0-9]+$ ]] && SPECIES="hs"
    [[ ${GENOME} =~ ^mm[0-9]+$ ]] && SPECIES="mm"
    [[ -z "${SPECIES}" ]] && echo "Unknown species for macs: ${GENOME}" && exit 1
    echo "${SPECIES}"
}