#!/usr/bin/env bash

echo "ChIP-Seq pipeline script"
WORK_DIR=`pwd`
echo "Working directory: $WORK_DIR"

# Load technical stuff
source ~/work/washu/scripts/util.sh

# Nothing to compare with aligned on hg38, so stick with hg19 for now
GENOME=hg19

INDEXES=${WORK_DIR}/../${GENOME}
~/work/washu/scripts/genome_indices.sh ${GENOME} ${INDEXES}
cd ${WORK_DIR}



echo "Submitting fastqc tasks"
FASTQC_TASKS=""
for FILE in $(find . -type f -name '*.f*q' -printf '%P\n')
do :
    QSUB_ID=`~/work/washu/scripts/fastqc.sh ${FILE}`
    echo "$FILE: $QSUB_ID"
    FASTQC_TASKS="$FASTQC_TASKS $QSUB_ID"
done
wait_complete ${FASTQC_TASKS}
check_logs
mkdir ${WORK_DIR}/fastqc
mv *_fastqc.* ${WORK_DIR}/fastqc
multiqc ${WORK_DIR}/fastqc



echo "Submitting trim 5 nucleotides tasks"
TRIM_TASKS=""
for FILE in $(find . -type f -name '*.f*q' -printf '%P\n')
do :
    QSUB_ID=`~/work/washu/scripts/trim.sh ${FILE} 5`
    echo "$FILE: $QSUB_ID"
    TRIM_TASKS="$TRIM_TASKS $QSUB_ID"
done
wait_complete ${TRIM_TASKS}
check_logs
TRIMMED=${WORK_DIR}/../trim
mkdir ${TRIMMED}
mv *_5.* ${TRIMMED}
cd ${TRIMMED}
WORK_DIR=`pwd`
echo "Working directory: $WORK_DIR"



echo "Submitting fastqc tasks"
FASTQC_TASKS=""
for FILE in $(find . -type f -name '*.f*q' -printf '%P\n')
do :
    QSUB_ID=`~/work/washu/scripts/fastqc.sh ${FILE}`
    echo "$FILE: $QSUB_ID"
    FASTQC_TASKS="$FASTQC_TASKS $QSUB_ID"
done
wait_complete ${FASTQC_TASKS}
check_logs
mkdir ${WORK_DIR}/fastqc
mv *_fastqc.* ${WORK_DIR}/fastqc
multiqc ${WORK_DIR}/fastqc



echo "Submitting bowtie tasks"
BOWTIE_TASKS=""
for FILE in $(find . -type f -name '*.f*q' -printf '%P\n')
do :
    QSUB_ID=`~/work/washu/scripts/bowtie.sh ${GENOME} ${FILE} ${INDEXES}`
    echo "$FILE: $QSUB_ID"
    BOWTIE_TASKS="$BOWTIE_TASKS $QSUB_ID"
done
wait_complete ${BOWTIE_TASKS}
check_logs
BAMS=${WORK_DIR}/../bams
mkdir ${BAMS}
mv *.bam ${BAMS}
mv *bowtie*.log ${BAMS}
cd ${BAMS}
WORK_DIR=`pwd`
echo "Working directory: $WORK_DIR"


echo "Create TDF tracks for $GENOME"
TDF_TASKS=""
for FILE in $(find . -type f -name '*.bam' -printf '%P\n')
do :
QSUB_ID=`~/work/washu/scripts/tdf.sh ${FILE} ${GENOME}`
    echo "$FILE: $QSUB_ID"
    SUBSAMPLE_TASKS="$TDF_TASKS $QSUB_ID"
done
wait_complete ${TDF_TASKS}
check_logs
mkdir ${WORK_DIR}/../tdfs
mv *tdf* ${WORK_DIR}/../tdfs


READS=15000000
echo "Subsampling to $READS reads"
SUBSAMPLE_TASKS=""
for FILE in $(find . -type f -name "*.bam" -printf '%P\n')
do :
    QSUB_ID=`~/work/washu/scripts/subsample.sh ${FILE} ${READS}`
    echo "$FILE: $QSUB_ID"
    SUBSAMPLE_TASKS="$SUBSAMPLE_TASKS $QSUB_ID"
done
wait_complete ${SUBSAMPLE_TASKS}
check_logs
SUBSAMPLED=${WORK_DIR}/../subsampled
mkdir ${SUBSAMPLED}
mv *${READS}* ${SUBSAMPLED}
cd ${SUBSAMPLED}
WORK_DIR=`pwd`
echo "Working directory: $WORK_DIR"



echo "Submitting macs2 tasks"
MACS2_TASKS=""
for FILE in $(find . -type f -name "*.bam" -printf '%P\n')
do :
    QSUB_ID=`~/work/washu/scripts/macs2.sh ${GENOME} 0.01 ${FILE}`
    echo "$FILE: $QSUB_ID"
    MACS2_TASKS="$MACS2_TASKS $QSUB_ID"
done
wait_complete ${MACS2_TASKS}
check_logs
PEAKS=${WORK_DIR}/../peaks
mkdir ${PEAKS}
mv *.bed ${PEAKS}
mv *macs* ${PEAKS}
cd ${PEAKS}
WORK_DIR=`pwd`
echo "Working directory: $WORK_DIR"