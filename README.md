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
https://genome.ucsc.edu/cgi-bin/hgFileUi?g=wgEncodeBroadHistone -> http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1003536

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

# ChIP-seq guidelines

ChIP-seq guidelines and practices of the ENCODE…
http://www.ncbi.nlm.nih.gov/pubmed/22955991

Sequencing and library complexity
For each ChIP-seq point-source library, ENCODE’s goal is to obtain 10 million uniquely mapping reads per replicate experiment for mammalian genomes, with a target NRF (nonredundancy fraction) 0.8 for 10 million reads.

IDR!

Practical Guidelines for the Comprehensive Analysis of ChIP-seq Data
https://scholar.google.com/scholar?hl=ru&q=Practical+Guidelines+for+the+Comprehensive+Analysis+of+ChIP-seq+Data&btnG=

# Reads QC

FastQC
http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Sequencing quality assessment tools to enable data-driven informatics for high throughput genomics
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3865868/#B1

NGS QC Toolkit: A Toolkit for Quality Control of Next Generation Sequencing Data
Quick comparison of existing tools for QC!
http://journals.plos.org/plosone/article/asset?id=10.1371%2Fjournal.pone.0030619.PDF

# Alignment

ENCODE: Mapping and analysis of chromatin state dynamics in nine human cell types
http://www.nature.com/nature/journal/v473/n7345/abs/nature09906.html
ChIP-seq reads were aligned to human genome build HG18 with MAQ (http://maq.sourceforge.net/maq-man.shtml) using default parameters. All reads were truncated to 36 bases before alignment.

Roadmapepigenomics: Integrative analysis of 111 reference human epigenomes
http://www.nature.com/nature/journal/v518/n7539/full/nature14248.html#methods
Sequenced data sets from the Release 9 of the Epigenome Atlas involved mapping a total of 150.21 billion sequencing reads onto hg19 assembly of the human genome using the PASH read mapper34. These read mappings were used (except for RNA-seq data sets, which were mapped as described above) for constructing the 111 consolidated epigenomes. Only uniquely mapping reads were retained and multiply-mapping reads were filtered out.
(Sampling to 30mln reads, ignoring black regions).

Major tools:
* bwa
* bowtie2
* pash
* MAQ

https://www.biostars.org/p/97197/
"I think it's important to note that BWA is one of the few fast mapping algorithms that allows for indels. Tools like Maq and Bowtie will not map reads if there is an insertion or deletion. I have used BWA to map 75bp Illumina reads at 20x coverage to a 30Mb fungal genome with good results."

# Peak calling
Practical Guidelines for the Comprehensive Analysis of ChIP-seq Data
http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003326

Major tools:
* MACS2
* SICER

# Differential Peak calling

A comprehensive comparison of tools for differential ChIP-seq analysis
http://bib.oxfordjournals.org/content/early/2016/01/12/bib.bbv110.abstract



# Pipeline
HTSeq—a Python framework to work with high-throughput sequencing data
https://bioinformatics.oxfordjournals.org/content/early/2014/10/20/bioinformatics.btu638.full
