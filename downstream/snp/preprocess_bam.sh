#!/usr/bin/env bash

if [ $# -lt 6 ]; then
    echo "Need 6 parameters! <bam> <out> <GID> <GPU> <GLB> <GSM>"
    exit 1
fi

READS_BAM=$1
OUTPUT_BAM=$2
GID=$3
GPU=$4
GLB=$5
GSM=$6

WORKING_DIR=$(dirname "$OUTPUT_BAM")

echo WORKING_DIR: ${WORKING_DIR}

cd ${WORKING_DIR}

FILE_NAME=$(basename "$READS_BAM")

NAME=${FILE_NAME%%.bam}

if [ ! -d "tmp" ]; then
	mkdir tmp
fi

echo Processing $NAME

samtools view -b -F 4 ${READS_BAM} >tmp/${NAME}_mapped.bam 

java -jar /mnt/stripe/tools/picard-tools/picard.jar AddOrReplaceReadGroups \
	I=tmp/${NAME}_mapped.bam \
	O=tmp/${NAME}_tmp.bam \
	RGID=${GID} \
	RGPU=${GPU} \
	RGSM=${GSM} \
	RGLB=${GLB} \
	RGPL=illumina || exit 1

samtools view -h tmp/${NAME}_tmp.bam | \
	grep  -v "chrUn_gl" | grep  -Ev "chr.+_gl.+random" | \
	grep  -Ev "chr17_ctg5_hap1|chr4_ctg9_hap1|chr6_apd_hap1|chr6_cox_hap2|chr6_dbb_hap3" | \
	grep  -Ev "chr6_mann_hap4|chr6_mcf_hap5|chr6_qbl_hap6|chr6_ssto_hap7" | \
	awk -v OFS='\t' '{ if ($5 == 255) $5=60; print}' | \
	samtools view -b >${OUTPUT_BAM} || exit 1

samtools index ${OUTPUT_BAM} || exit 1

rm tmp/${NAME}_mapped.bam
rm tmp/${NAME}_tmp.bam
