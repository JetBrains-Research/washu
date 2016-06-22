#!/usr/bin/env bash

echo "ChIP-Seq pipeline script"
WORK_DIR=`pwd`
echo "Working directory: $WORK_DIR"

# Load technical stuff
source ~/work/washu/scripts/util.sh

# Nothing to compare with aligned on hg38, so stick with hg19 for now
GENOME=hg19

INDEXES=${WORK_DIR}/../${GENOME}
echo "Genomes and indices folder: ${INDEXES}"
~/work/washu/scripts/genome_indices.sh ${GENOME} ${INDEXES}
cd ${WORK_DIR}

# Batch QC
~/work/washu/scripts/fastqc.sh ${WORK_DIR}

# Batch trim
~/work/washu/scripts/trim.sh ${WORK_DIR} 5

# Move trimmed out to _trim folder and CD
TRIMMED=${WORK_DIR}_trim
mkdir ${TRIMMED}
mv *_5.* ${TRIMMED}
cd ${TRIMMED}
WORK_DIR=`pwd`
echo "Working directory: $WORK_DIR"

# Batch QC in trimmed folder
~/work/washu/scripts/fastqc.sh ${WORK_DIR}

# Batch Bowtie
~/work/washu/scripts/bowtie.sh ${WORK_DIR} ${GENOME} ${INDEXES}

# Move results and CD
BAMS=${WORK_DIR}_bams
mkdir ${BAMS}
mv *.bam ${BAMS}
mv *bowtie*.log ${BAMS}
cd ${BAMS}
WORK_DIR=`pwd`
echo "Working directory: $WORK_DIR"


# Batch TDF
~/work/washu/scripts/tdf.sh ${WORK_DIR} ${GENOME}

# Move results
mkdir ${WORK_DIR}_tdfs
mv *.tdf ${WORK_DIR}_tdfs
mv *tdf.log ${WORK_DIR}_tdfs


# Batch subsampling
READS=15000000
~/work/washu/scripts/subsample.sh ${WORK_DIR} ${READS}

# Move results and CD
SUBSAMPLED=${WORK_DIR}_subsampled
mkdir ${SUBSAMPLED}
mv *${READS}* ${SUBSAMPLED}
cd ${SUBSAMPLED}
WORK_DIR=`pwd`
echo "Working directory: $WORK_DIR"


# Batch macs2
~/work/washu/scripts/macs2.sh ${WORK_DIR} ${GENOME} 0.01

# Move results
PEAKS=${WORK_DIR}_peaks
mkdir ${PEAKS}
mv *.bed ${PEAKS}
mv *macs* ${PEAKS}

# Batch zinbra
~/work/washu/scripts/zinbra.sh ${WORK_DIR} ${GENOME} ${INDEXES} 0.01

# Move results
ZINBRA_PEAKS=${WORK_DIR}_zinbra_peaks
mkdir ${ZINBRA_PEAKS}
mv *.bed ${ZINBRA_PEAKS}

echo "Done"