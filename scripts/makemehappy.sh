#!/usr/bin/env bash

echo "Pipeline script1"
echo "Working directory: `pwd`"

echo "Downloading files"
~/work/washu/scripts/download.sh

echo "Submitting sra to fastq.gz tasks"
find . -type f -name "*.sra" | xargs -n1 ~/work/washu/scripts/sra2fastq.sh

# See sra2fastq script for tasks naming convention
echo "Collecting tasks: sra2fastq"
SRA_TASKS=""
for FILE in $(find . -type f -name "*.sra")
do :
    NAME=${FILE%%.sra} # file name without extension
    if [ ! -n "$SRA_TASKS" ]; then
        SRA_TASKS="sra2fastq_$NAME"
    else
        SRA_TASKS="$SRA_TASKS,sra2fastq_$NAME"
    fi
done

echo "Waiting for tasks $SRA_TASKS to finish"
qsub -hold_jid $SRA_TASKS -cwd ~/work/washu/scripts/makemehappy2.sh