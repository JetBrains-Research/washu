#!/usr/bin/env bash

# Check tools
which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
which macs2 &>/dev/null || { echo "macs2 not found! Install macs2: <https://github.com/taoliu/MACS/wiki/Install-macs2>"; exit 1; }

# Load technical stuff
source $(dirname $0)/../parallel/util.sh
SCRIPT_DIR="$(project_root_dir)"
export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p "${TMPDIR}"

################################################################################
# Configuration start ##########################################################
################################################################################

>&2 echo "run_chipdiff: $@"
if [ $# -lt 4 ]; then
    echo "Need 4 parameters! <NAME> <WORKING_DIR> <CHROM_SIZES> <DIFF_MACS_POOLED>"
    exit 1
fi

NAME=$1

CHIPDIFF=$2

CHROM_SIZES=$3

DIFF_MACS_POOLED=$4

BROAD_CUTOFF=0.1
echo "BROAD_CUTOFF $BROAD_CUTOFF"

READS_Y=$(find . -name 'YD*.bam' | sed 's#\./##g' | grep -v 'input' | tr '\n' ' ')

READS_O=$(find . -name 'OD*.bam' | sed 's#\./##g' | grep -v 'input' | tr '\n' ' ')

# Check MACS2 for shift values
macs2_shift() {
    echo $(cat $1 | grep "# d =" | sed 's/.*# d = //g')
}

echo
echo "Processing $CHIPDIFF"
echo $CHIPDIFF

echo "Processing chipdiff as on pooled tags (reads)"

cat >config.txt <<CONFIG
maxIterationNum  500
minP             0.95
maxTrainingSeqNum 10000
minFoldChange    3
minRegionDist    1000
CONFIG

>&2 echo "Processing Y Tags"
SHIFT_Y=$(macs2_shift ${DIFF_MACS_POOLED}/Y_${BROAD_CUTOFF}_peaks.xls)
echo "SHIFT Y: $SHIFT_Y"

run_parallel << SCRIPT
#!/bin/sh
#PBS -N Y_bam2tags
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=8gb
#PBS -j oe
#PBS -o ${CHIPDIFF}/Y_bams2tags.log
# This is necessary because qsub default working dir is user home
cd ${CHIPDIFF}
module load bedtools2
for F in ${READS_Y}; do
    >&2 echo \$F
    bash ${SCRIPT_DIR}/scripts/bam2tags.sh \$F $SHIFT_Y >> Y_tags.tag
done
SCRIPT

QSUB_ID1=$QSUB_ID

>&2 echo "Processing O Tags"
SHIFT_O=$(macs2_shift ${DIFF_MACS_POOLED}/O_${BROAD_CUTOFF}_peaks.xls)
echo "SHIFT O: $SHIFT_O"
run_parallel << SCRIPT
#!/bin/sh
#PBS -N O_bam2tags
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=8gb
#PBS -j oe
#PBS -o ${CHIPDIFF}/O_bams2tags.log
# This is necessary because qsub default working dir is user home
cd ${CHIPDIFF}
module load bedtools2
for F in ${READS_O}; do
    >&2 echo \$F
    bash ${SCRIPT_DIR}/scripts/bam2tags.sh \$F $SHIFT_O >> O_tags.tag
done
SCRIPT
    QSUB_ID2=$QSUB_ID

wait_complete "$QSUB_ID1 $QSUB_ID2"

run_parallel << SCRIPT
#!/bin/sh
#PBS -N ${NAME}_chipdiff_3
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${CHIPDIFF}/${NAME}_chipdiff_3.log
# This is necessary because qsub default working dir is user home
cd ${CHIPDIFF}

source "${SCRIPT_DIR}/parallel/util.sh"
export TMPDIR=\$(type job_tmp_dir &>/dev/null && echo "\$(job_tmp_dir)" || echo "/tmp")

# Inplace sort
sort -T \${TMPDIR} -k1,1 -k2,2n -o Y_tags.tag Y_tags.tag
sort -T \${TMPDIR} -k1,1 -k2,2n -o O_tags.tag O_tags.tag

ChIPDiff Y_tags.tag O_tags.tag $CHROM_SIZES config.txt ${NAME}_3
cat ${NAME}_3.region | awk -v OFS='\t' '\$4=="-" {print \$1,\$2,\$3}' > ${NAME}_3_cond1.bed
cat ${NAME}_3.region | awk -v OFS='\t' '\$4=="+" {print \$1,\$2,\$3}' > ${NAME}_3_cond2.bed
SCRIPT

wait_complete $QSUB_ID

rm Y_tags.tag O_tags.tag

# TMP dir cleanup:
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir
