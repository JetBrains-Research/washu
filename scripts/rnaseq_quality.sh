#!/usr/bin/env bash
# original rna-seq quality script by Alex Predeus
# modified by zayats1812@mail.ru

# Load technical stuff
source ~/work/washu/scripts/util.sh

WORK_DIR=$1

RRNA="/scratch/artyomov_lab_aging/indexes/hg19/hg19.rRNA_merged.intervals"
REFFLAT="/scratch/artyomov_lab_aging/indexes/hg19/hg19.refFlat.txt"
GENOME="/scratch/artyomov_lab_aging/indexes/hg19/GRCh37.p13.genome.fa"

PICARD="/home/kzaytsev/epiProject/tools/picardTools/picard-tools-2.4.1/picard.jar"


echo "Batch rna-seq quality: ${WORK_DIR}"

cd ${WORK_DIR}

TASKS=""
for FILE in $(find . -type f -name '*.bam' -printf '%P\n' | grep -vE ".tr.")
do :
    TAG=${FILE%%.bam} # file name without extension

    # Submit task
    QSUB_ID=$(qsub -v WORK_DIR="$WORK_DIR",TAG="$TAG" ~/work/washu/scripts/rnaseq_quality_task.sh)
    echo "FILE: ${FILE}; JOB: ${QSUB_ID}"
    TASKS="$TASKS $QSUB_ID"
done
wait_complete ${TASKS}

echo -e "Sample\tN_reads\tPct_mapped\tPct_mapped_1loc\tPct_unmapped\tPct_rRNA\tPct_coding\tPct_UTR\tPct_intronic\tPct_intergenic\tJunctions\tInsertion_rate\tDeletion_rate\tPct_NC_junctions\tDel_av_length\tIns_av_length" > allSamples.rnastat
cat ./*.rnastat >> allSamples.rnastat

echo "Done. Batch rna-seq quality: ${WORK_DIR}"