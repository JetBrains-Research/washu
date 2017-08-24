#!/usr/bin/env bash

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

export PATH=/opt/conda/envs/py3.5/bin:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/SICER_V1.1/SICER

# Redirect both stderr and stdout to stdout then tee and then to stderr
bash $QSUB_FILE &> "$LOG"