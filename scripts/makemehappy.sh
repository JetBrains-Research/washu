#!/usr/bin/env bash

# Small routing to wait until all the tasks are finished on the qsub cluster
wait_complete()
{
    JOBS=$1
    TITLE=$2
    echo "Waiting: $TITLE $JOBS"
    for JOB in $JOBS
    do :
        # The job id is actually the first numbers in the string
        JOB_ID=`echo $JOB | awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}'`
        if [ ! -z "$JOB_ID" ]; then
            while qstat $JOB_ID &> /dev/null; do
                echo -n "."
                sleep 5
            done;
        fi
        echo
    done
    echo "Finished: $TITLE"
}

echo "Pipeline script"
echo "Working directory: `pwd`"

echo "Downloading files"
~/work/washu/scripts/download.sh

echo "Submitting sra2fastq tasks"
SRA2FASTQ_JOBS=""
for FILE in $(find . -type f -name "*.sra")
do :
    QSUB_ID=`~/work/washu/scripts/sra2fastq.sh $FILE`
    SRA2FASTQ_JOBS="$SRA2FASTQ_JOBS $QSUB_ID"
done
wait_complete $SRA2FASTQ_JOBS "sra2fastq tasks"

echo "Submitting bowtie tasks"
BOWTIE_JOBS=""
for FILE in $(find . -type f -name "*.sra")
do :
    QSUB_ID=`~/work/washu/scripts/bowtie.sh hg37 $FILE`
    BOWTIE_JOBS="$BOWTIE_JOBS $QSUB_ID"
done
wait_complete $BOWTIE_JOBS "bowtie tasks"