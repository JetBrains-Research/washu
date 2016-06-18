#!/usr/bin/env bash

echo "ChIP-Seq pipeline script"
WORK_DIR=`pwd`
echo "Working directory: $WORK_DIR"

echo "Load util stuff"

# Small procedure to wait until all the tasks are finished on the qsub cluster
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

# Checks for logs
check_logs()
{

    ERRORS=`find . -name "*.log" | xargs grep -i -e "err|warn"`
    if [ ! -z "$ERRORS" ]; then
        echo "ERRORS found"
        echo "$ERRORS"
        exit 1
    fi
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
    chmod a+r *
    cd $WORK_DIR
fi

echo "Submitting sra2fastq if tasks if necessary"
SRA2FASTQ_JOBS=""
for FILE in $(find . -type f -name "*.sra" -printf '%P\n')
do :
    QSUB_ID=`~/work/washu/scripts/sra2fastq.sh $FILE`
    echo "$FILE: $QSUB_ID"
    SRA2FASTQ_JOBS="$SRA2FASTQ_JOBS $QSUB_ID"
done
wait_complete $SRA2FASTQ_JOBS
check_logs

echo "Submitting fastqc tasks"
FASTQC_TASKS=""
for FILE in $(find . -type f -regextype posix-egrep -regex '.*(fastq|fq)(\.gz)?' -printf '%P\n')
do :
    QSUB_ID=`~/work/washu/scripts/fastqc.sh $FILE`
    echo "$FILE: $QSUB_ID"
    FASTQC_TASKS="$FASTQC_TASKS $QSUB_ID"
done
wait_complete $FASTQC_TASKS
check_logs

echo "Check bowtie indexes"
if [ ! -f "$WORK_DIR/$GENOME/$GENOME.1.ebwt" ]; then
    cd $WORK_DIR/$GENOME
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N bowtie_indexes_${GENOME}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o $WORK_DIR/qsub/bowtie_indexes_${GENOME}.log

# Load module
module load bowtie

# This is necessary because qsub default working dir is user home
cd $WORK_DIR/$GENOME
bowtie-build $(find $WORK_DIR/$GENOME -type f -name "*.fa" -printf '%P\n' | paste -sd "," -) $GENOME
ENDINPUT
)
    cd $WORK_DIR
    wait_complete $QSUB_ID
    check_logs
fi

echo "Submitting bowtie tasks"
BOWTIE_TASKS=""
for FILE in $(find . -type f -regextype posix-egrep -regex '.*(fastq|fq)' -printf '%P\n')
do :
    QSUB_ID=`~/work/washu/scripts/bowtie.sh $GENOME $FILE`
    echo "$FILE: $QSUB_ID"
    BOWTIE_TASKS="$BOWTIE_TASKS $QSUB_ID"
done
wait_complete $BOWTIE_TASKS
check_logs

READS=15000000
echo "Subsampling to $READS reads"
SUBSAMPLE_TASKS=""
for FILE in $(find . -type f -name "*.bam" -printf '%P\n')
do :
    QSUB_ID=`~/work/washu/scripts/subsample.sh $FILE $READS`
    echo "$FILE: $QSUB_ID"
    SUBSAMPLE_TASKS="$SUBSAMPLE_TASKS $QSUB_ID"
done
wait_complete $SUBSAMPLE_TASKS
check_logs

echo "Submitting macs2 tasks"
MACS2_TASKS=""
for FILE in $(find . -type f -name "*$READS*.bam" -printf '%P\n')
do :
    QSUB_ID=`~/work/washu/scripts/macs2.sh $GENOME 0.01 $FILE`
    echo "$FILE: $QSUB_ID"
    MACS2_TASKS="$MACS2_TASKS $QSUB_ID"
done
wait_complete $MACS2_TASKS
check_logs