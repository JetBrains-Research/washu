License [![license](https://img.shields.io/github/license/mashape/apistatus.svg)](https://opensource.org/licenses/MIT)
Tests [![tests](http://teamcity.jetbrains.com/app/rest/builds/buildType:(id:Epigenome_Tools_Washu)/statusIcon)](http://teamcity.jetbrains.com/viewType.html?buildTypeId=Epigenome_Tools_Washu&guest=1)
Pipeline tests [![long tests](http://teamcity.jetbrains.com/app/rest/builds/buildType:(id:Epigenome_Tools_WashuPipelineTests)/statusIcon)](http://teamcity.jetbrains.com/viewType.html?buildTypeId=Epigenome_Tools_WashuPipelineTests&guest=1)  

Technical pipelines
===================
Technical pipelines for ChIP-Seq, RNA-Seq, RRBS processing on Portable Batch System (qsub).

Project
-------

* `/bed`            - BED files manipulations - intersection, ChromHMM enrichment, etc.
* `/docker`         - Docker configuration files with tools and test data. Used for `pipeline_chipseq.py` testing.
* `/notebooks`      - Jupiter notebooks
* `/parallel`       - Scripts for parallel execution of PBS greed using `qsub` queue management
* `/R`              - R scripts
* `/reports`        - Scripts to generate report by logs (bowtie, macs2, etc) 
* `/scripts`        - Visualization, downstream analysis, e.g. replicated ChIP-Seq comparison, etc.
* `/test`           - Tests
* `/tex`            - Latex presentation for the `pipeline_chipseq.py`
* `/uscs`           - Prepared custom tracks for UCSC genome [browser](https://genome.ucsc.edu/)
* `/web`            - Create web server with tracks and peaks powered by [Biodalliance](http://www.biodalliance.org/) genome browser.


Scripts
-------
Scripts for batch processing.
* `fastqc.sh`   - Quality control for raw-reads using [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 
and [MultiQC](http://multiqc.info/) for aggregated report  
* `bowtie.sh`   - Bowtie alignment of all the read files in a folder. Creates summary report:

| sample |  reads | aligned | not_aligned | supressed |
| ------------- |-------------:| -----:| -----:| -----:|
| CD14_H3K4me3_hg38_bowtie.log| 50000000 | 49000000  (98.0%) | 300000 (0.6%) | 7000000 (1.4%) |  
* `bigwig.sh`   - Alignment files visualization suitable for UCSC in a folder
* `macs2.sh`    - Peak calling using [MACS2](https://github.com/taoliu/MACS). 
Creates summary report, with (**R**eads **I**n **P**eaks) and (**F**raction of **R**eads **I**n **P**eaks) metrics

| sample |  tags | redundant_rate | paired_peaks | fragment| alternatives | peaks| rip| frip |
| ------------- |-------------:| -----:| -----:| -----:|-----:|-----:|-----:|-----:|
| CD14_H3K4me3_hg38_bowtie_macs_broad_0.1.log| 50000000 | 0.1 | 20500 | 150 | 3,150 | 40000 | 3500000 | 70 |
* `rseg.sh`     - Peak calling using [RSeg](https://academic.oup.com/bioinformatics/article/27/6/870/236489/Identifying-dispersed-epigenomic-domains-from-ChIP) algorithm
* `sicer.sh`    - Broad histone marks peak calling using [SICER](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2732366/) algorithm
* `rpkm.sh`     - Create deduplicated normalized tracks using [RPKM](http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/)				

All the scripts are designed for batch processing, and print all the necessary tools and arguments.

Parallel execution
------------------
All the scripts from `/paralles` are suitable to work with Portable Batch System and on single machine.
Parallelism level on local machine can be configured via **WASHU_PARALLELISM** environment variable.
 
Pipelines
---------
* `pipeline_chipseq.py`         - Pipeline for batch ULI-ChIP-Seq processing, including QC, alignment, peak calling
* `pipeline_utils.py`           - Pipeline for batch RNA-Seq data processing
* `scripts/chipseq_diff.sh`     - Pipeline for replicated ChIP-Seq comparison using MACS2, [DiffBind](http://www.nature.com/nature/journal/v481/n7381/full/nature10730.html), 
[ChIPDiff](https://academic.oup.com/bioinformatics/article/24/20/2344/258202/An-HMM-approach-to-genome-wide-identification-of), 
[MANorm](https://www.ncbi.nlm.nih.gov/pubmed/22424423)


Requirements
------------
Add the following to `~/.bashrc` (Linux) or `~/.bash_profile` (MacOS):
```bash
# Allow pipeline execution from anywhere
export PYTHONPATH="<PATH_TO_REPOSITORY>:$PYTHONPATH"
```

Docker
------
There is a Docker hub image `biolabs/washu` with all the necessary data and tools for pipeline and tests.
Docker configurations are available under `/docker` folder.

Update necessary docker image.
```bash
docker pull biolabs/washu
```


Tests
-----
```bash
# FAST tests
docker run -v <project_path>:/washu -t -i biolabs/washu /bin/bash -c "source activate py3.5 && cd /washu && bash test.sh"

# SLOW tests 
docker run -v <project_path>:/washu -t -i biolabs/washu /bin/bash -c "source activate py3.5 && cd /washu && bash test_pipeline.sh"
```

Data standards and pipelines
--------------
* ENCODE data [standards](https://www.encodeproject.org/data-standards/)
* ENCODE ChIP-Seq [pipeline](https://github.com/ENCODE-DCC/chip-seq-pipeline)
* Blueprint ChIP-Seq [pipeline](http://dcc.blueprint-epigenome.eu/#/md/chip_seq_grch38)

Useful links
------------
* Washington University in Saint Louis Maxim Artyomov LAB [homepage](https://artyomovlab.wustl.edu/site/) 
* JetBrains BioLabs [homepage](https://research.jetbrains.org/groups/biolabs)
* Review on ChIP-Seq, ATAC-Seq and DNAse-Seq processing in latex [format](https://github.com/olegs/bioinformatics)