#!/usr/bin/env bash

echo "Pipeline script"
echo "Working directory: `pwd`"

echo "Load util stuff"
# Small routing to wait until all the tasks are finished on the qsub cluster
wait_complete()
{
    echo "Waiting for tasks."
    for JOB in $@
    do :
        echo -n "JOB: $JOB"
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
    echo "Done."
}

echo "Downloading files"
~/work/washu/scripts/download.sh

echo "Submitting sra2fastq tasks"
SRA2FASTQ_JOBS=""
for FILE in $(find . -type f -name "*.sra")
do :
    QSUB_ID=`~/work/washu/scripts/sra2fastq.sh $FILE`
    SRA2FASTQ_JOBS="$SRA2FASTQ_JOBS $QSUB_ID"
done
wait_complete $SRA2FASTQ_JOBS

echo "Submitting bowtie tasks"
BOWTIE_TASKS=""
for FILE in $(find . -type f -name "*.fastq")
do :
    QSUB_ID=`~/work/washu/scripts/bowtie.sh hg19 $FILE`
    BOWTIE_TASKS="$BOWTIE_TASKS $QSUB_ID"
done
wait_complete $BAM2BED_TASKS

echo "Submitting bam2bed tasks"
BAM2BED_TASKS=""
for FILE in $(find . -type f -name "*.bam")
do :
    QSUB_ID=`~/work/washu/scripts/bam2bed.sh $FILE`
    BAM2BED_TASKS="$BAM2BED_TASKS $QSUB_ID"
done
wait_complete $BAM2BED_TASKS