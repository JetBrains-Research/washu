#!/bin/bash
# author oleg.shpynov@jetbrains.com
# Submit multiple tasks with xargs:
# find . -name "*.sra" -print0 | xargs -r0 -n1 sra2fastq.sh

WORK_DIR=`pwd`
SRA_FILE=$1

echo "Converting sra to fastq-dump $SRA_FILE"

qsub << ENDINPUT
#!/bin/sh
#PBS -N sra2fastq_$SRA_FILE
#PBS -l nodes=1:ppn=1,walltime=2:00:00,vmem=6gb
#PBS -j oe
#PBS -q dque

# Loading sratoolkit module
module load sratoolkit

cd $WORK_DIR
fastq-dump --split-3 --outdir $WORK_DIR $SRA_FILE
ENDINPUT
