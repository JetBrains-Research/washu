#!/usr/bin/env bash
# author oleg.shpynov@jetbrains.com

which STAR &>/dev/null || { echo "STAR not found! Download STAR: <https://github.com/alexdobin/STAR>"; exit 1; }

# Check configuration
[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util/util.sh

>&2 echo "index-star $@"
if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <GENOME> <FOLDER>"
    exit 1
fi
GENOME=$1
FOLDER=$2

STAR_INDEX_FOLDER="${FOLDER}/star"
echo "STAR index folder: ${STAR_INDEX_FOLDER}"
if [ -f "${STAR_INDEX_FOLDER}/SAindex" ]; then
    echo "Indexes already exist at ${STAR_INDEX_FOLDER}"
    exit 0
fi

GFT_FILE="${STAR_INDEX_FOLDER}"/${GENOME}.gtf

mkdir -p "${STAR_INDEX_FOLDER}"
cd "${STAR_INDEX_FOLDER}"
GTF_URL=""
if [ ${GENOME} = "mm9" ]; then
    GTF_URL="ftp://anonymous@ftp.ensembl.org/pub/release-67/gtf/mus_musculus/Mus_musculus.NCBIM37.67.gtf.gz"
elif [ ${GENOME} = "mm10" ]; then
    GTF_URL="ftp://anonymous@ftp.ensembl.org/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.gtf.gz"
elif [ ${GENOME} = "hg18" ]; then
    GTF_URL="ftp://anonymous@ftp.ensembl.org/pub/release-54/gtf/homo_sapiens/Homo_sapiens.NCBI36.54.gtf.gz"
elif [ ${GENOME} = "hg19" ]; then
    GTF_URL="ftp://anonymous@ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz"
elif [ ${GENOME} = "hg38" ]; then
    GTF_URL="ftp://anonymous@ftp.ensembl.org/pub/release-89/gtf/homo_sapiens/Homo_sapiens.GRCh38.89.gtf.gz"
else
    echo "Unknown genome ${GENOME} build, failed to download GTF file"
    exit 1
fi
echo "Downloading GTF file for build ${GENOME} ${GTF_URL}"
wget -nc "${GTF_URL}"
gunzip *.gz
DOWNLOADED_GTF=$(ls ${STAR_INDEX_FOLDER}/*.gtf)

# Add chr prefix if necessary
CHR_PREFIX=$(head -n 10 ${DOWNLOADED_GTF} | grep '^chr')
if [ -z "${CHR_PREFIX}" ]; then
    grep '^#!ge' ${DOWNLOADED_GTF} > ${GFT_FILE}
    grep -v '^#!ge' ${DOWNLOADED_GTF} | awk '{printf("chr%s\n", $0)}' >> ${GFT_FILE}
else
    cp ${DOWNLOADED_GTF} ${GFT_FILE}
fi
echo "GTF_FILE = ${GFT_FILE}"



if [ ! -f "${GFT_FILE}" ]; then
    echo "Failed to find GTF file for ${GENOME} in ${STAR_INDEX_FOLDER}"
    exit 1;
fi


cd ${FOLDER}
run_parallel << SCRIPT
#!/bin/sh
#PBS -N star_indexes_${GENOME}
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=64gb
#PBS -j oe
#PBS -o ${FOLDER}/${GENOME}_star_indexes.log

# This is necessary because qsub default working dir is user home
cd ${STAR_INDEX_FOLDER}

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ${STAR_INDEX_FOLDER} \
    --genomeFastaFiles ${FOLDER}/*.fa \
    --sjdbGTFfile ${GFT_FILE} \
    --sjdbOverhang 100

SCRIPT

wait_complete ${QSUB_ID}
check_logs

>&2 echo "Done. index-star $@"