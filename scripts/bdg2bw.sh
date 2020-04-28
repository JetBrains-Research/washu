#!/bin/bash

# this script is from Tao Liu https://gist.github.com/taoliu/2469050 

which bedtools &>/dev/null || { echo "ERROR: bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/> or conda install bedtools"; exit 1; }
which bedGraphToBigWig &>/dev/null || { echo "ERROR: bedGraphToBigWig not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/> or conda install -c bioconda ucsc-bedgraphtobigwig"; exit 1; }
which bedClip &>/dev/null || { echo "ERROR: bedClip not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/> or conda install -c bioconda ucsc-bedclip"; exit 1; }
 
# end of checking
 
if [ $# -lt 2 ]; then
    echo "Need 2 parameters! <bedgraph> <chrom info> [output]"
    exit 1
fi
 
BDG_FILE=$1
CHROM_SIZES=$2
BW_FILE=$3

# Check configuration
[[ ! -z ${WASHU_ROOT} ]] || { echo "ERROR: WASHU_ROOT not configured"; exit 1; }
source ${WASHU_ROOT}/parallel/util.sh
export TMPDIR=$(type job_tmp_dir &>/dev/null && echo "$(job_tmp_dir)" || echo "/tmp")
mkdir -p "${TMPDIR}"
if [ -z "$BW_FILE" ]; then
  BW_FILE=${BDG_FILE/.bdg/.bw}
fi

echo "Converting ${BDG_FILE} -> ${BW_FILE}"
# Remove coordinates outside chromosome sizes
bedtools slop -i ${BDG_FILE} -g ${CHROM_SIZES} -b 0 | bedClip stdin ${CHROM_SIZES} ${BDG_FILE}.clip
# Fix problem with not sorted clip file
LC_COLLATE=C sort -k1,1 -k2,2n -T ${TMPDIR} ${BDG_FILE}.clip > ${BDG_FILE}.sort.clip
bedGraphToBigWig ${BDG_FILE}.sort.clip ${CHROM_SIZES} ${BW_FILE}

# Cleanup
rm -f ${BDG_FILE}*.clip
type clean_job_tmp_dir &>/dev/null && clean_job_tmp_dir
