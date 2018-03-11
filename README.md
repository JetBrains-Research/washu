License [![license](https://img.shields.io/github/license/mashape/apistatus.svg)](https://opensource.org/licenses/MIT)
Tests [![tests](http://teamcity.jetbrains.com/app/rest/builds/buildType:(id:Epigenome_Tools_Washu)/statusIcon)](http://teamcity.jetbrains.com/viewType.html?buildTypeId=Epigenome_Tools_Washu&guest=1)
Pipeline tests [![long tests](http://teamcity.jetbrains.com/app/rest/builds/buildType:(id:Epigenome_Tools_WashuPipelineTests)/statusIcon)](http://teamcity.jetbrains.com/viewType.html?buildTypeId=Epigenome_Tools_WashuPipelineTests&guest=1)  

Technical pipelines
===================
Technical pipelines for ChIP-Seq, RNA-Seq, RRBS processing on Portable Batch System (qsub).

How do I launch the pipeline?
--------------------------
Follow these instructions to launch ChIP-Seq pipeline:
* Configure environment, see **Requirements** section
* Place all the `.fastq` files to a single `<FASTQ_FOLDER>`
* Create `<INDEXES>` folder to store all the indexes required
* Launch the pipeline with desired `<genome>`, e.g. `mm9` or `hg19` 
```bash
python pipeline_chipseq.py <FASTQ_FOLDER> <INDEXES> <genome>
```

Requirements
------------
* Ensure you have Python 3 installed as default interpreter
* Add the following to `~/.bashrc` (Linux) or `~/.bash_profile` (MacOS):
```bash
# Configure project path
export WASHU_ROOT="<PATH_TO_REPOSITORY>"

# Configure correct python code execution
export PYTHONPATH="$WASHU_ROOT:$PYTHONPATH"

# Configure local machine parallelism
export WASHU_PARALLELISM=8
```

* Install required tools using [Conda](https://conda.io/docs/)
```bash
conda install --channel bioconda samtools bedtools bowtie bowtie2 fastqc sra-tools macs2 \
    deeptools star rseg   ucsc-bedgraphtobigwig ucsc-bedclip ucsc-bigwigaverageoverbed
```
* Download [Picard tools](https://github.com/broadinstitute/picard):
```bash 
curl --location https://github.com/broadinstitute/picard/releases/download/2.10.7/picard.jar \
    --output ~/picard.jar
```
* Download and extract [Phantompeakqualtools](https://github.com/kundajelab/phantompeakqualtools):
```bash
curl --location https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/phantompeakqualtools/ccQualityControl.v.1.1.tar.gz \
    --output ~/phantompeakqualtools.tar.gz 
tar xvf ~/phantompeakqualtools.tar.gz
```
Project
-------

* `/bed`            - BED files manipulations - intersection, ChromHMM enrichment, closes gene, etc.
* `/docker`         - Docker configuration files with tools and test data. See Tests section.
* `/downstream`     - Aging project downstream analysis
* `/parallel`       - Scripts for parallel execution of PBS greed using `qsub` queue management
* `/scripts`        - QC, Visualization, BAM conversions, Reads In Peaks, etc.
* `/test`           - Tests

Scripts
-------
Scripts for batch processing.
* `fastqc.sh`   - Quality control for raw-reads using [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 
and [MultiQC](http://multiqc.info/) for aggregated report  
* `bowtie.sh`   - [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml) alignment of all the read files in a folder. Creates summary report:

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
* `downstream`                  - Downstream analysis including differential ChIP-Seq analysis

Docker
------
Build Docker image `biolabs/washu` with all the necessary data and tools for pipeline and tests.
Docker configurations are available under `/docker` folder.

Ensure your `biolabs/test-data` image is up to date.
```bash
docker pull biolabs/test-data
```

Build image:
```bash
docker build -t biolabs/washu /docker/biolabs/washu
```

Tests
-----
```bash
# Tests
docker run -v <project_path>:/washu -t -i biolabs/washu /bin/bash -c "source activate py3.5 && cd /washu && bash test.sh"

# Pipeline tests 
docker run -v <project_path>:/washu -t -i biolabs/washu /bin/bash -c "source activate py3.5 && cd /washu && bash test_pipeline.sh"
```

Tools used
---------- 
* [Bedtools](https://bedtools.readthedocs.io/en/latest/)
* [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml)
* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [ChIPDiff](https://academic.oup.com/bioinformatics/article/24/20/2344/258202/An-HMM-approach-to-genome-wide-identification-of) 
* [Deeptools](http://deeptools.readthedocs.io/en/latest/content/installation.html)
* [DiffBind](http://www.nature.com/nature/journal/v481/n7381/full/nature10730.html) 
* [DiffReps](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0065598)
* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 
* [MACS1.4](https://github.com/taoliu/MACS)
* [MACS2](https://github.com/taoliu/MACS)
* [MANorm](https://www.ncbi.nlm.nih.gov/pubmed/22424423)
* [MultiQC](http://multiqc.info/)
* [Phantompeakqualtools](https://github.com/kundajelab/phantompeakqualtools)
* [Picardtools](https://github.com/broadinstitute/picard)
* [RSeg](https://academic.oup.com/bioinformatics/article/27/6/870/236489/Identifying-dispersed-epigenomic-domains-from-ChIP)
* [RSem](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323)
* [Samtools](http://samtools.sourceforge.net/)
* [SICER](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2732366/)
* [Star](https://www.ncbi.nlm.nih.gov/pubmed/23104886)
* [ZINBRA](https://github.com/JetBrains-Research/zinbra)

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