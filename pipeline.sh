#!/usr/bin/env bash

# Configuration
WORK_DIR=`pwd`
GENOME=hg19 # Nothing to compare with aligned on hg38
INDEXES=${WORK_DIR}/../${GENOME}

# Load technical stuff
source ~/work/washu/scripts/util.sh

echo "Genomes and indices folder: ${INDEXES}"
bash ~/work/washu/scripts/genome_indices.sh ${GENOME} ${INDEXES}
cd ${WORK_DIR}

# Batch QC
bash ~/work/washu/scripts/fastqc.sh ${WORK_DIR}

# Batch trim
bash ~/work/washu/scripts/trim.sh ${WORK_DIR} 5

# Move trimmed out to _trim folder and CD
TRIMMED=${WORK_DIR}_trim
mkdir ${TRIMMED}
mv *_5.* ${TRIMMED}
cd ${TRIMMED}
WORK_DIR=`pwd`
echo "Working directory: $WORK_DIR"

# Batch QC in trimmed folder
bash ~/work/washu/scripts/fastqc.sh ${WORK_DIR}

# Batch Bowtie
bash ~/work/washu/scripts/bowtie.sh ${WORK_DIR} ${GENOME} ${INDEXES}

# Move results and CD
BAMS=${WORK_DIR}_bams
mkdir ${BAMS}
mv *.bam ${BAMS}
mv *bowtie*.log ${BAMS}
cd ${BAMS}
# Create summary
python ~/work/washu/bowtie_logs.py ${BAMS}
WORK_DIR=`pwd`
echo "Working directory: $WORK_DIR"


# Batch TDF
bash ~/work/washu/scripts/tdf.sh ${WORK_DIR} ${GENOME}

# Move results
mkdir ${WORK_DIR}_tdfs
mv *.tdf ${WORK_DIR}_tdfs
mv *tdf.log ${WORK_DIR}_tdfs


#####################
# Batch subsampling #
#####################
READS=15000000
bash ~/work/washu/scripts/subsample.sh ${WORK_DIR} ${READS}

# Move results and CD
SUBSAMPLED=${WORK_DIR}_subsampled
mkdir ${SUBSAMPLED}
mv *${READS}* ${SUBSAMPLED}
cd ${SUBSAMPLED}
WORK_DIR=`pwd`
echo "Working directory: $WORK_DIR"


# Batch macs with different peak calling procedures settings
QS=( 0.001 0.01 0.1)
for Q in "${QS[@]}"
do
    bash ~/work/washu/scripts/macs2.sh ${WORK_DIR} ${GENOME} ${Q}
    # Move results
    PEAKS=${WORK_DIR}_macs_${Q}
    mkdir ${PEAKS}
    mv *.bed ${PEAKS}
    mv *macs* ${PEAKS}
    # Create summary
    python ~/work/washu/macs2_logs.py ${PEAKS}
done

echo "Done"