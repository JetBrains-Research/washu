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
    if [[ ! ${FILE} =~ ^.*input.*$ ]]; then
        # WashU data naming convention
        >&2 echo "NOT INPUT ${FILE}"
        DONOR=$(echo ${FILE} | sed -e "s/.*\(donor[^_\.]*\).*/\1/")
        if [[ ! -z ${DONOR} ]]; then
            >&2 echo "DONOR ${DONOR}"
            INPUTS=$(find . -name "*${DONOR}.*input*.bam" -printf '%P\n')
            INPUT=${INPUTS[0]}
        fi

        # ENCODE UW or Broad
        if [[ ${FILE} =~ ^Broad.*_1_.*$ ]]; then
            >&2 echo "Broad file _1"
            INPUTS=$(find . -name "Broad*_input_1*.bam" -printf '%P\n')
            INPUT=${INPUTS[0]}
        fi
        if [[ ${FILE} =~ ^Broad.*_2_.*$ ]]; then
            >&2 echo "Broad file _1"
            INPUTS=$(find . -name "Broad*_input_2*.bam" -printf '%P\n')
            INPUT=${INPUTS[0]}
        fi
        if [[ ${FILE} =~ ^UW.+$ ]]; then
            >&2 echo "UW file"
            INPUTS=$(find . -name "UW*_input*.bam" -printf '%P\n')
            INPUT=${INPUTS[0]}
        fi
        >&2 echo "FILE ${FILE} INPUT ${INPUT}"
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