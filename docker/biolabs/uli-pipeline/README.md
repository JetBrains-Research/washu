Docker Image with ULI-ChIP-seq pipeline tools
=====================================

This is just an image `ubuntu:latest` - Ubuntu LTS with all the environment setup and necessary tools.

Build Docker image
---------

With `washu` repository cloned, navigate to `<washu>/docker/biolabs/uli-pipeline` and run `docker build -t biolabls/uli-pipeline .`.

```bash
<washu>$ cd ./docker/biolabs/uli-pipeline
<washu>/docker/biolabs/uli-pipeline$ docker build -t biolabs/uli-pipeline .
```

Building the image will take some time.


Prepare data
---------

Choose a working directory.
* Make sure the directory is writable.
* Make sure there's enough disk space free (the necessary amount depends on how much data you intend to process).

Make subdirectories named `fastq` and `indexes`.

```bash
<working_directory>$ mkdir fastq indexes
```

Place your FASTQ files in `fastq` subdirectory.
* If you have input files, rename them so that their name contains `input`.

```bash
<working_directory>/fastq$ mv SRR1234567890.fastq SRR1234567890_input.fastq
```

If you have prebuilt bowtie indexes for your genome build, place them in `indexes/<build>`.
* If you don't have the indexes, they will be generated automatically.

```bash
<working_directory>/indexes$ mkdir mm9
<working_directory>/indexes$ cp ~/bowtie/mm9* ./mm9/
```

Launch pipeline
----------

Note the location of your working directory and of the cloned `washu` repository.
* The working directory is the parent directory of `fastq` and `indexes` subfolders.

Run the following command:
```bash
docker run -v <washu>:/washu -v <working_directory>:/data -e WASHU_ROOT=/washu -e LOCAL_USER_ID=`id -u $USER` -it biolabs/uli-pipeline /bin/bash -c "bash /washu/pipeline_chipseq.sh <genome_build>"
```

The pipeline should take from several hours to several days depending on the available computing capabilities and the amount of data.

Example:
```bash
$ docker run -v /home/user/washu:/washu -v /home/user/work:/data -e WASHU_ROOT=/washu -e LOCAL_USER_ID=`id -u $USER` -it biolabs/uli-pipeline /bin/bash -c "bash /washu/pipeline_chipseq.sh mm9"
```
