#!/usr/bin/env bash

echo "Pipeline script"
echo "Working directory: `pwd`"

echo "Downloading files"
~/work/washu/scripts/download.sh

echo "Submitting sra to fastq.gz tasks"
find . -type f -name "*.sra" | xargs -n1 ~/work/washu/scripts/sra2fastq.sh

echo "Waiting for sra2fastq tasks to finish"
for FILE in $(find . -type f -name "*.sra")
do :
    NAME=${FILE%%.sra} # file name without extension
    JOB_ID=$(qsub "sra2fastq_$NAME") # See sra2fastq script for tasks naming convention
    while qstat $JOB_ID &> /dev/null; do
        echo -n "."
        sleep 5
    done;
    echo
done
echo "All sra2fastq tasks are finished"