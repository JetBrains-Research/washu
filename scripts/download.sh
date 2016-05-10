#!/bin/bash
# Download H3K4me3 data for Monocytes CD14+ from ENCODE and roadmapepigenomics
# author oleg.shpynov@jetbrains.com

echo "Downloading ENCODE data"
if [ ! -f SRR568364.sra ]; then
    wget -r ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR568/SRR568364/SRR568364.sra -O SRR568364.sra
fi
if [ ! -f SRR568365.sra ]; then
    wget -r ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR568/SRR568365/SRR568365.sra -O SRR568365.sra
fi
echo "Downloading Roadmapepigenomics data"
if [ ! -f SRR787515.sra ]; then
    wget -r ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR787/SRR787515/SRR787515.sra -O SRR787515.sra
fi
if [ ! -f SRR787516.sra ]; then
    wget -r ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR787/SRR787516/SRR787516.sra -O SRR787516.sra
fi