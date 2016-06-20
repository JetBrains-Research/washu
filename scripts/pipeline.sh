#!/usr/bin/env bash

echo "ChIP-Seq pipeline script"
WORK_DIR=`pwd`
echo "Working directory: $WORK_DIR"

# Load technical stuff
source util.sh

echo "Log folder $WORK_DIR/qsub"
if [ ! -d "$WORK_DIR/qsub" ]; then
    mkdir $WORK_DIR/qsub
fi

# Check indices
GENOME=hg38
./genome_indices.sh $GENOME $WORK_DIR
BOWTIE_INDEXES=${WORK_DIR}/${GENOME}



echo "Submitting fastqc tasks"
FASTQC_TASKS=""
for FILE in $(find . -type f -name '*.f*q' -printf '%P\n')
do :
    QSUB_ID=`~/work/washu/scripts/fastqc.sh $FILE`
    echo "$FILE: $QSUB_ID"
    FASTQC_TASKS="$FASTQC_TASKS $QSUB_ID"
done
wait_complete $FASTQC_TASKS
check_logs
mkdir ${WORK_DIR}/fastqc
mv *_fastqc.* ${WORK_DIR}/fastqc



echo "Submitting trim 5 nucleotides tasks"
TRIM_TASKS=""
for FILE in $(find . -type f -name '*.f*q' -printf '%P\n')
do :
    QSUB_ID=`~/work/washu/scripts/trim.sh $FILE 5`
    echo "$FILE: $QSUB_ID"
    TRIM_TASKS="$TRIM_TASKS $QSUB_ID"
done
wait_complete $TRIM_TASKS
check_logs
mkdir ${WORK_DIR}/../trim
mv *_5.fq ${WORK_DIR}/../trim
cd ${WORK_DIR}/../trim
WORK_DIR=`pwd`
echo "Working directory: $WORK_DIR"



echo "Submitting fastqc tasks"
FASTQC_TASKS=""
for FILE in $(find . -type f -name '*.f*q' -printf '%P\n')
do :
    QSUB_ID=`~/work/washu/scripts/fastqc.sh $FILE`
    echo "$FILE: $QSUB_ID"
    FASTQC_TASKS="$FASTQC_TASKS $QSUB_ID"
done
wait_complete $FASTQC_TASKS
check_logs
mkdir ${WORK_DIR}/fastqc
mv *_fastqc.* ${WORK_DIR}/fastqc



echo "Submitting bowtie tasks"
BOWTIE_TASKS=""
for FILE in $(find . -type f -name '.f*q' -printf '%P\n')
do :
    QSUB_ID=`~/work/washu/scripts/bowtie.sh $GENOME $FILE $BOWTIE_INDEXES`
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
mkdir ${WORK_DIR}/../subsampled
mv *${READS}* ${WORK_DIR}/../subsampled
cd ${WORK_DIR}/../subsampled
WORK_DIR=`pwd`
echo "Working directory: $WORK_DIR"



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