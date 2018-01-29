#!/usr/bin/env bash
# original rna-seq quality script by Alex Predeus
# modified by zayats1812@mail.ru
# TODO: fix hardcoded!

# Check configuration
[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util/util.sh

>&2 echo "Batch rnaseq-quality $@"
if [ $# -lt 1 ]; then
    echo "Need 1 parameter! <WORK_DIR>"
    exit 1
fi
WORK_DIR=$1

RRNA="/scratch/artyomov_lab_aging/indexes/hg19/hg19.rRNA_merged.intervals"
REFFLAT="/scratch/artyomov_lab_aging/indexes/hg19/hg19.refFlat.txt"
GENOME="/scratch/artyomov_lab_aging/indexes/hg19/GRCh37.p13.genome.fa"

cd ${WORK_DIR}

TASKS=()
for FILE in $(find . -type f -name '*.bam' | sed 's#\./##g' | grep -vE ".tr.")
do :
    TAG=${FILE%%.bam} # file name without extension

    # Submit task
    QSUB_ID=$(qsub -v WORK_DIR="$WORK_DIR",TAG="$TAG" ${WASHU_ROOT}/parallel/rnaseq_quality_task.sh)
    echo "FILE: ${FILE}; TASK: ${QSUB_ID}"
    TASKS+=("$QSUB_ID")
done
wait_complete ${TASKS[@]}

echo -e "Sample\tN_reads\tPct_mapped\tPct_mapped_1loc\tPct_unmapped\tPct_rRNA\tPct_coding\tPct_UTR\tPct_intronic\tPct_intergenic\tJunctions\tInsertion_rate\tDeletion_rate\tPct_NC_junctions\tDel_av_length\tIns_av_length" > allSamples.rnastat
cat ./*.rnastat >> allSamples.rnastat

echo "Done. Batch rna-seq quality: ${WORK_DIR}"
