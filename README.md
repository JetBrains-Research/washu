Technical pipelines for ChIP-Seq, RNA-Seq, RRBS processing on PBS greed
==========

Code
----
* `/R`          - R scripts
* `/bed`        - Scripts for bed files manipulation - intersection, chromhmm enrichment, etc.
* `/logs`       - Scripts for batch logs (bowtie, macs2, etc) processing 
* `/notebooks`  - Jupiter notebooks
* `/scripts`    - Various scripts for execution of PBS greed using `qsub` queue management
* `/test`       - Tests
* `/uscs`       - Prepared custom tracks for UCSC genome [browser](https://genome.ucsc.edu/)

Pipelines
---------

* `pipeline_chipseq.py` - Pipeline for batch ULI-ChIP-Seq processing, including QC, alignment, peak calling
* `pipeline_utils.py`   - Pipeline for batch RNA-Seq data processing


Useful links
------------
* Washington University in Saint Louis Maxim Artyomov LAB [homepage](https://artyomovlab.wustl.edu/site/) 
* JetBrains BioLabs [homepage](https://research.jetbrains.org/groups/biolabs)