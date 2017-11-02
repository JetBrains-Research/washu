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

>&2 echo "chipseq_diff: $@"
if [ $# -lt 3 ]; then
    echo "Need 3 parameters! <NAME> <CHROM_SIZES> <DIFFBIND_CSV>"
    exit 1
fi

NAME=$1
CHROM_SIZES=$2
DIFFBIND_CSV=$3

FOLDER=$(pwd)
echo "FOLDER"
echo $FOLDER

PREFIX="$FOLDER/$NAME"
echo "PREFIX"
echo $PREFIX

Q=0.01
echo "Q $Q"
BROAD_CUTOFF=0.1
echo "BROAD_CUTOFF $BROAD_CUTOFF"

READS_Y=$(awk -v FS=',' '{ if ($4 == "Y") print $6 }' $DIFFBIND_CSV | sort  -T ${TMPDIR} --unique | tr '\n' ' ')
echo "READS Y"
echo "$READS_Y"
READS_O=$(awk -v FS=',' '{ if ($4 == "O") print $6 }' $DIFFBIND_CSV | sort -T ${TMPDIR} --unique | tr '\n' ' ')
echo "READS O"
echo "$READS_O"

INPUTS_Y=$(awk -v FS=',' '{ if ($4 == "Y") print $8 }' $DIFFBIND_CSV | sort -T ${TMPDIR} --unique | tr '\n' ' ')
echo "INPUT_READS Y"
echo "$INPUTS_Y"
INPUTS_O=$(awk -v FS=',' '{ if ($4 == "O") print $8 }' $DIFFBIND_CSV | sort -T ${TMPDIR} --unique | tr '\n' ' ')
echo "INPUT_READS O"
echo "$INPUTS_O"

################################################################################
# Configuration end ############################################################
################################################################################

DIFF_MACS_POOLED="${PREFIX}_macs_pooled"
echo
echo "Processing $DIFF_MACS_POOLED"
if [ ! -d $DIFF_MACS_POOLED ]; then
    mkdir -p ${DIFF_MACS_POOLED}
    cd ${DIFF_MACS_POOLED}
    echo "Processing MACS2 pooled peaks and compare them";

    run_parallel << SCRIPT
#!/bin/sh
#PBS -N ${NAME}_Y_macs2_broad_${BROAD_CUTOFF}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${DIFF_MACS_POOLED}/${NAME}_Y_macs2_broad_${BROAD_CUTOFF}.log
# This is necessary because qsub default working dir is user home
cd ${DIFF_MACS_POOLED}
macs2 callpeak -t $READS_Y -c $INPUTS_Y -f BAM -g hs -n Y_${BROAD_CUTOFF} -B --broad --broad-cutoff ${BROAD_CUTOFF}
SCRIPT
    QSUB_ID1=$QSUB_ID

    run_parallel << SCRIPT
#!/bin/sh
#PBS -N ${NAME}_O_macs2_broad_${BROAD_CUTOFF}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${DIFF_MACS_POOLED}/${NAME}_O_macs2_broad_${BROAD_CUTOFF}.log
# This is necessary because qsub default working dir is user home
cd ${DIFF_MACS_POOLED}
macs2 callpeak -t $READS_O -c $INPUTS_O -f BAM -g hs -n O_${BROAD_CUTOFF} -B --broad --broad-cutoff ${BROAD_CUTOFF}
SCRIPT
    QSUB_ID2=$QSUB_ID
    wait_complete "$QSUB_ID1 $QSUB_ID2"

    check_logs
    bash ${SCRIPT_DIR}/bed/compare.sh Y_${BROAD_CUTOFF}_peaks.broadPeak O_${BROAD_CUTOFF}_peaks.broadPeak ${NAME}_${BROAD_CUTOFF}
fi


# MANorm
MANORM="${PREFIX}_manorm"
echo
echo "Processing $MANORM"
if [ ! -d $MANORM ]; then
    mkdir -p ${MANORM}

    echo "Processing MAnorm using pooled MACS2 peaks as peakfile and pooled reads as readfiles"
# README.txt
# Create a folder and place in the folder MAnorm.sh, MAnorm.r, and all 4 bed files to be analyzed.
# run command:   ./MAnorm.sh    sample1_peakfile[BED]     sample2_peakfile[BED] \
#                               sample1_readfile[BED]     sample2_readfile[BED]  \
#                               sample1_readshift_lentgh[INT]      sample2_readshift_length[INT]
    MANORM_SH=$(which MAnorm.sh)
    echo "Found MAnorm.sh: ${MANORM_SH}"
    cp ${MANORM_SH} ${MANORM_SH%%.sh}.r ${MANORM}

    cp ${DIFF_MACS_POOLED}/Y_${BROAD_CUTOFF}_peaks.broadPeak ${MANORM}/Y_peaks.bed
    cp ${DIFF_MACS_POOLED}/O_${BROAD_CUTOFF}_peaks.broadPeak ${MANORM}/O_peaks.bed

    >&2 echo "Processing Y Pooled Reads"
    run_parallel << SCRIPT
#!/bin/sh
#PBS -N Y_bam2reads
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=8gb
#PBS -j oe
#PBS -o ${MANORM}/Y_bams2reads.log
# This is necessary because qsub default working dir is user home
cd ${MANORM}

module load bedtools2
for F in ${READS_Y}; do
    >&2 echo \$F
    bash ${SCRIPT_DIR}/scripts/bam2reads.sh \$F >> Y_reads.bed
done
SCRIPT
    QSUB_ID1=$QSUB_ID

    >&2 echo "Processing O Pooled Reads"
    run_parallel << SCRIPT
#!/bin/sh
#PBS -N O_bam2reads
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=8gb
#PBS -j oe
#PBS -o ${MANORM}/O_bams2reads.log
# This is necessary because qsub default working dir is user home
cd ${MANORM}

module load bedtools2
for F in ${READS_O}; do
    >&2 echo \$F
    bash ${SCRIPT_DIR}/scripts/bam2reads.sh \$F >> O_reads.bed
done
SCRIPT
    QSUB_ID2=$QSUB_ID

    wait_complete "$QSUB_ID1 $QSUB_ID2"

    # Check MACS2 for shift values
    macs2_shift() {
        echo $(cat $1 | grep "# d =" | sed 's/.*# d = //g')
    }

    SHIFT_Y=$(macs2_shift ${DIFF_MACS_POOLED}/Y_${BROAD_CUTOFF}_peaks.xls)
    echo "SHIFT Y: $SHIFT_Y"

    SHIFT_O=$(macs2_shift ${DIFF_MACS_POOLED}/O_${BROAD_CUTOFF}_peaks.xls)
    echo "SHIFT O: $SHIFT_O"

    run_parallel << SCRIPT
#!/bin/sh
#PBS -N manorm_${Q}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=16gb
#PBS -j oe
#PBS -o ${MANORM}/manorm_${Q}.log
# This is necessary because qsub default working dir is user home
cd ${MANORM}

source "${SCRIPT_DIR}/parallel/util.sh"
export TMPDIR=\$(type job_tmp_dir &>/dev/null && echo "\$(job_tmp_dir)" || echo \"/tmp\")

# Sort inplace
sort -T \${TMPDIR} -k1,1 -k2,2n -o Y_reads.bed Y_reads.bed
sort -T \${TMPDIR} -k1,1 -k2,2n -o Y_peaks.bed Y_peaks.bed

sort -T \${TMPDIR} -k1,1 -k2,2n -o O_reads.bed O_reads.bed
sort -T \${TMPDIR} -k1,1 -k2,2n -o O_peaks.bed O_peaks.bed

# Load required modules
module load R
module load bedtools2

bash ${MANORM}/MAnorm.sh Y_peaks.bed O_peaks.bed Y_reads.bed O_reads.bed $SHIFT_Y $SHIFT_O
SCRIPT
    wait_complete $QSUB_ID
fi

# TMP dir cleanup:
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir
