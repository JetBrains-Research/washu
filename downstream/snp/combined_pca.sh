#!/usr/bin/env bash


if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <SNP_FILE> <G1000_FILE>"
    exit 1
fi

GATK_TOOL="/mnt/stripe/tools/gatk-4.0.3.0/gatk"
PLINK_TOOL="/mnt/stripe/tools/PLINK/plink"
SNP_FILE=$1
g1000=$2

gunzip <${SNP_FILE} | awk '{gsub(/^chr/,""); print}' | awk '{gsub(/ID=chr/,"ID="); print}' | bgzip >no_chr.vcf.gz

tabix -p vcf no_chr.vcf.gz

${GATK_TOOL} SelectVariants \
        -select "vc.isBiallelic() && AF > 0.05" \
        -V ${g1000} \
        -O s1000.vcf.gz || exit 1

${GATK_TOOL} SelectVariants \
        --concordance no_chr.vcf.gz \
        -V s1000.vcf.gz \
        -O conc_1000.vcf.gz || exit 1

bcftools merge no_chr.vcf.gz conc_1000.vcf.gz | bgzip >merged.vcf.gz

tabix -p vcf merged.vcf.gz || exit 1

${GATK_TOOL} SelectVariants \
        --max-nocall-number 0 \
        -V merged.vcf.gz \
        -O final.vcf.gz || exit 1

${PLINK_TOOL} --allow-extra-chr --double-id --vcf final.vcf.gz --out output

${PLINK_TOOL} --allow-extra-chr --bfile output --pca 39
