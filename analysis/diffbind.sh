#!/usr/bin/env bash
# Script to process differential chip-seq analysis
# author Oleg Shpynov (oleg.shpynov@jetbrains.com)

# Check tools
which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
which R &>/dev/null || { echo "R not found! Install R <https://www.r-project.org/>"; exit 1; }
>&2 echo "diffbind: $@"

if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <DIFFBIND_CSV> <GENES_GTF>"
    exit 1
fi

CSV=$1
CSV_NAME=${CSV%%.csv}
GENES_GTF=$2

# Load cluster stuff
source ~/work/washu/scripts/util.sh

################################################################################
# Configuration start ##########################################################
################################################################################

GROUP1="Y"
echo "GROUP1"
echo $GROUP1

GROUP2="O"
echo "GROUP2"
echo $GROUP2

NAME="diff_k27ac_${GROUP1}_${GROUP2}"
echo "NAME"
echo $NAME

echo "Use diffbind"
if [ ! -f "${CSV_NAME}_result.csv" ]; then
    WORK_DIR=$(pwd)
    QSUB_ID=$(qsub << ENDINPUT
#!/bin/sh
#PBS -N ${NAME}_diffbind
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${NAME}_diffbind.log
# This is necessary because qsub default working dir is user home
cd ${WORK_DIR}
module load R
Rscript ~/work/washu/analysis/diffbind.R ${CSV}
ENDINPUT
)
    wait_complete "$QSUB_ID"
    check_logs
fi

# Filter out old and young donors and sort by Q-Value
cat ${CSV_NAME}_result.csv | awk 'NR > 1 {print $0}' | sed 's#"##g' | tr ',' '\t' |\
 awk -v OFS='\t' '$9 > 0 {print $0}' | sort -k1,1 -k2,2n > ${CSV_NAME}_cond1.bed
cat ${CSV_NAME}_result.csv | awk 'NR > 1 {print $0}' | sed 's#"##g' | tr ',' '\t' |\
 awk -v OFS='\t' '$9 < 0 {print $0}' | sort -k1,1 -k2,2n > ${CSV_NAME}_cond2.bed

# Save diffbind results to simple BED3 format
awk -v OFS='\t' '{ print $1,$2,$3}' ${CSV_NAME}_cond1.bed > ${CSV_NAME}_cond1.bed3
awk -v OFS='\t' '{ print $1,$2,$3}' ${CSV_NAME}_cond2.bed > ${CSV_NAME}_cond2.bed3

~/work/washu/bed/closest_gene.sh ${GENES_GTF} \
    ${CSV_NAME}_cond1.bed3 > ${CSV_NAME}_cond1_closest_genes.tsv
~/work/washu/bed/closest_gene.sh ${GENES_GTF} \
    ${CSV_NAME}_cond2.bed3 > ${CSV_NAME}_cond2_closest_genes.tsv
