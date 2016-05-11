#!/usr/bin/env bash

echo "Pipeline script"
WORK_DIR=`pwd`
echo "Working directory: $WORK_DIR"

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

echo "Log folder $WORK_DIR/qsub"
if [ ! -d "$WORK_DIR/qsub" ]; then
    mkdir $WORK_DIR/qsub
fi

GENOME=hg38
echo "Processing genome sequence $GENOME"
if [ ! -d "$WORK_DIR/$GENOME" ]; then
    mkdir $WORK_DIR/$GENOME
fi
if [ ! -f "$WORK_DIR/$GENOME/chr1.fa" ]; then
    cd $WORK_DIR/$GENOME
    # Download only chromosomes sequences
    rsync -avzP --exclude="chr*_*" --exclude="*.txt" rsync://hgdownload.cse.ucsc.edu/goldenPath/$GENOME/chromosomes/ .
    gunzip *.fa.gz
    cd $WORK_DIR
fi

echo "Downloading files"
~/work/washu/scripts/download.sh

echo "Submitting sra2fastq tasks"
SRA2FASTQ_JOBS=""
for FILE in $(find . -type f -name "*.sra" -printf '%P\n')
do :
    QSUB_ID=`~/work/washu/scripts/sra2fastq.sh $FILE`
    SRA2FASTQ_JOBS="$SRA2FASTQ_JOBS $QSUB_ID"
done
wait_complete $SRA2FASTQ_JOBS

echo "Submitting fastqc tasks"
FASTQC_TASKS=""
for FILE in $(find . -type f -name "*.fastq" -printf '%P\n')
do :
    QSUB_ID=`~/work/washu/scripts/fastqc.sh $FILE`
    FASTQC_TASKS="$FASTQC_TASKS $QSUB_ID"
done
wait_complete $FASTQC_TASKS

echo "Check bowtie indexes"
if [ ! -f "$WORK_DIR/$GENOME/$GENOME.1.ebwt" ]; then
    cd $WORK_DIR/$GENOME
    # Load module
    module load bowtie
    bowtie-build $(find . -type f -name "*.fa" -printf '%P\n' | paste -sd "," -) $GENOME
    cd $WORK_DIR
fi

echo "Submitting bowtie tasks"
BOWTIE_TASKS=""
for FILE in $(find . -type f -name "*.fastq" -printf '%P\n')
do :
    QSUB_ID=`~/work/washu/scripts/bowtie.sh $GENOME $FILE`
    BOWTIE_TASKS="$BOWTIE_TASKS $QSUB_ID"
done
wait_complete $BOWTIE_TASKS

echo "Submitting macs2 tasks"
MACS2_TASKS=""
for FILE in $(find . -type f -name "*.bam" -printf '%P\n')
do :
    QSUB_ID=`~/work/washu/scripts/macs2.sh $GENOME 0.01 $FILE`
    MACS2_TASKS="$MACS2_TASKS $QSUB_ID"
done
wait_complete $MACS2_TASKS