#!/usr/bin/env bash
# This is a chip-seq technical pipeline scratch.
#
# Contains the following steps:
#   * FastQC
#   * Trim 5
#   * FastQC
#   * Alignment with BWA, summary for alignment
#   * Subsampling to 15mln reads
#   * BigWigs visualization
#   * Peak calling MACS2 narrow peaks, fraction of reads in peaks, summary
#   * Peak calling MACS2 broad peaks, fraction of reads in peaks, summary
#
# Usage:
#   * Launch FastQC on the whole bunch of fq files.
#   * Decide whether trim or subsampling is required
#   * Modify inplace copy of this pipeline. (cp ~/work/washu/pipeline.sh <working_folder_with_fq>)
#   * Launch pipeline, and wait for "Done" message.
#
# Conventions:
# This pipeline uses folder naming as a steps. I.e. next step appends _suffix for the working folder,
# stores results in new folder and change working folder if necessary.
# All the indices are stored in <working_folder>/../<genome>.
#
# Example:
# run_6_7 -> run_6_7_trim -> run_6_7_trim_bams -> run_6_7_trim_bams_bws
#                                              -> run_6_7_trim_bams_macs_0.01
#                                              -> run_6_7_trim_bams_macs_broad_0.01
#
# author oleg.shpynov@jetbrains.com


# Configuration
WORK_DIR=`pwd`
GENOME=hg19 # Nothing to compare with aligned on hg38
INDEXES=${WORK_DIR}/../${GENOME}

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
# MultiQC is able to process Bowtie report
echo "Processing multiqc"
multiqc ${WORK_DIR}

# Batch BigWig visualization
bash ~/work/washu/scripts/bigwig.sh ${WORK_DIR} ${INDEXES}/${GENOME}.chrom.sizes
# Move results
mkdir ${WORK_DIR}_bws
mv *.bw ${WORK_DIR}_bws
mv *.bdg ${WORK_DIR}_bws
mv *bw.log ${WORK_DIR}_bws

# Batch subsampling
READS=15
echo "Subsampling to ${READS}mln"
bash ~/work/washu/scripts/subsample.sh ${WORK_DIR} ${READS}

# Move results and CD
SUBSAMPLED=${WORK_DIR}_${READS}mln
mkdir ${SUBSAMPLED}
mv *${READS}* ${SUBSAMPLED}
cd ${SUBSAMPLED}
WORK_DIR=`pwd`
echo "Working directory: $WORK_DIR"

# Batch BigWig visualization
bash ~/work/washu/scripts/bigwig.sh ${WORK_DIR} ${INDEXES}/${GENOME}.chrom.sizes
# Move results
mkdir ${WORK_DIR}_bws
mv *.bw ${WORK_DIR}_bws
mv *.bdg ${WORK_DIR}_bws
mv *bw.log ${WORK_DIR}_bws

# Batch macs with different peak calling procedures settings
QS=( 0.001 0.01 0.1)
for Q in "${QS[@]}"
do
    bash ~/work/washu/scripts/macs2.sh ${WORK_DIR} ${GENOME} ${Q} ${INDEXES}/${GENOME}.chrom.sizes
    # Move results
    PEAKS=${WORK_DIR}_macs_${Q}
    mkdir ${PEAKS}
    mv *${Q}* ${PEAKS}
    # Create summary
    python ~/work/washu/macs2_logs.py ${PEAKS}
done

# Batch macs with different peak calling procedures settings
QS=( 0.001 0.01 0.1)
for Q in "${QS[@]}"
do
    bash ~/work/washu/scripts/macs2_broad.sh ${WORK_DIR} ${GENOME} ${Q} ${INDEXES}/${GENOME}.chrom.sizes
    # Move results
    PEAKS=${WORK_DIR}_macs_broad_${Q}
    mkdir ${PEAKS}
    mv *${Q}* ${PEAKS}
    # Create summary
    python ~/work/washu/macs2_logs.py ${PEAKS}
done

echo "Done"