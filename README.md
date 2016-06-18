# ChIP-Seq Pipeline

Technical pipeline for ChIP-Seq processing on PBS greed. 
Most of the computational intensive steps are being executed in parallel on the greed.

### Input
* Genome build for analysis (hg38) - hardcoded.

### Steps
* Download genome sequence
* SRA -> fastq (if necessary)
* Quality control (fastqc)
* Build indexes for aligner (bowtie)
* Map reads on reference genome
* Call narrow peaks (macs2)

NOTE: this is the very first sketch, tools used are subject to discuss.

### Usage
* Navigate to directory shared among the greed
* Execute `/scripts/chipseq.sh` script

# Data
Please ensure that you have all the data downloaded.
