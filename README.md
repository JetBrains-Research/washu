Technical pipelines for ChIP-Seq, RNA-Seq, RRBS processing on Portable Batch System (qsub)
==========

Project
-------
* `/R`          - R scripts
* `/analysis`   - Downstream analysis, e.g. replicated ChIP-Seq comparison
* `/bed`        - Scripts for bed files manipulation - intersection, ChromHMM enrichment, etc.
* `/reports`    - Scripts to generate report by logs (bowtie, macs2, etc) 
* `/notebooks`  - Jupiter notebooks
* `/scripts`    - Various scripts for execution of PBS greed using `qsub` queue management
* `/test`       - Tests
* `/uscs`       - Prepared custom tracks for UCSC genome [browser](https://genome.ucsc.edu/)
* `/web`        - Create web server with tracks and peaks powered by [Biodalliance](http://www.biodalliance.org/) genome browser
Scripts
-------
Useful scripts for batch processing. Suitable to work with Portable Batch System and on single machine
* `fastqc.sh`   - Quality control for raw-reads using [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MultiQC](http://multiqc.info/) for aggregated report  
* `bowtie.sh`   - Bowtie alignment of all the read files in a folder. Creates summary report:

| sample |  reads | aligned | not_aligned | supressed |
| ------------- |-------------:| -----:| -----:| -----:|
| CD14_H3K4me3_hg38_bowtie.log| 50000000 | 49000000  (98.0%) | 300000 (0.6%) | 7000000 (1.4%) |  
* `bigwig.sh`   - Alignment files visualization suitable for UCSC in a folder
* `macs2.sh`    - Peak calling using [MACS2](https://github.com/taoliu/MACS). Creates summary report, with (**R**eads **I**n **P**eaks) and (**F**raction of **R**eads **I**n **P**eaks) metrics.
* `rseg.sh`     - Peak calling using [RSeg](https://academic.oup.com/bioinformatics/article/27/6/870/236489/Identifying-dispersed-epigenomic-domains-from-ChIP) algorithm
* `sicer.sh`    - Broad histone marks peak calling using [SICER](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2732366/) algorithm
| sample |  tags | redundant_rate | paired_peaks | fragment| alternatives | peaks| rip| frip |
| ------------- |-------------:| -----:| -----:| -----:|-----:|-----:|-----:|-----:|
| CD14_H3K4me3_hg38_bowtie_macs_broad_0.1.log| 50000000 | 0.1 | 20500 | 150 | 3,150 | 40000 | 3500000 | 70 |
* `rpkm.sh`     - Create deduplicated normalized tracks using [RPKM](http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/)				

All the scripts are designed for batch processing, and print all the necessary tools and arguments.

Pipelines
---------
* `pipeline_chipseq.py`         - Pipeline for batch ULI-ChIP-Seq processing, including QC, alignment, peak calling
* `pipeline_utils.py`           - Pipeline for batch RNA-Seq data processing
* `analysis/chipseq_diff.sh`    - Pipeline for replicated ChIP-Seq comparison using MACS2, [DiffBind](http://www.nature.com/nature/journal/v481/n7381/full/nature10730.html), 
[ChIPDiff](https://academic.oup.com/bioinformatics/article/24/20/2344/258202/An-HMM-approach-to-genome-wide-identification-of), 
[MANorm](https://www.ncbi.nlm.nih.gov/pubmed/22424423)


Data standards and pipelines
--------------
* ENCODE data [standards](https://www.encodeproject.org/data-standards/)
* ENCODE ChIP-Seq [pipeline](https://github.com/ENCODE-DCC/chip-seq-pipeline)
* Blueprint ChIP-Seq [pipeline](http://dcc.blueprint-epigenome.eu/#/md/chip_seq_grch38)

Useful links
------------
* Washington University in Saint Louis Maxim Artyomov LAB [homepage](https://artyomovlab.wustl.edu/site/) 
* JetBrains BioLabs [homepage](https://research.jetbrains.org/groups/biolabs)
* Review on ChIP-Seq, ATAC-Seq and DNAse-Seq processing in latex format: https://github.com/olegs/bioinformatics