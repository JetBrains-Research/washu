# washu
WashU project

# Pipeline

Technical pipeline for ChIP-Seq processing on PBS greed. 
Most of the computational intensive steps are being executed in parallel on the greed.

Input:
* Genome build for analysis (hg38) - hardcoded.

Steps:
* Download genome sequence
* Download data SRA reads
* SRA -> fastq
* Quality control (fastqc)
* Build indexes for aligner (bowtie)
* Map reads on reference genome
* Call narrow peaks (macs2)

NOTE: this is the very first sketch, tools used are subject to discuss.

Usage:
* Navigate to directory shared among the greed
* Execute `/scripts/makemehappy.sh` script
* Profit!

# Data
Monocytes CD14+; 2 donors; H3K4me3

# Data ENCODE
https://genome.ucsc.edu/cgi-bin/hgFileUi?g=wgEncodeBroadHistone -> 
http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1003536

Single strand reads
* Run1
http://sra.dnanexus.com/runs/SRR568364
* Run2
http://sra.dnanexus.com/runs/SRR568365


# Data Roadmap Epigenomics
http://www.ncbi.nlm.nih.gov/geo/roadmap/epigenomics/?display=500 ->
http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1102797

Single strand reads
* Run1
http://sra.dnanexus.com/runs/SRR787515
* Run2
http://sra.dnanexus.com/runs/SRR787516

