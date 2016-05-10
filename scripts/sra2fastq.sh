#!/bin/bash
# author oleg.shpynov@jetbrains.com

WORK_DIR=`pwd`
SRA_FILE=$1

echo "Converting sra to fastq-dump $SRA_FILE"

qsub << ENDINPUT
#!/bin/sh
#PBS -N sra2fastq$SRA_FILE
#PBS -l nodes=1:ppn=1,walltime=2:00:00,vmem=6gb
#PBS -j oe
#PBS -q dque

cd $WORK_DIR
fastq-dump --log-level err --dumpbase --gzip --outdir $WORK_DIR $SRA_FILE
ENDINPUT
