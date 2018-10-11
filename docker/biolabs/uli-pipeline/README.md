Docker Image with ULI-ChIP-seq pipeline tools
=====================================

This is just an image `ubuntu:latest` - Ubuntu LTS with all the environment setup and necessary tools.

Build Docker image
---------

With `washu` repository cloned, navigate to `<washu>/docker/biolabs/uli-pipeline` and run `sudo docker build -t biolabls/uli-pipeline .`.
You might need to provide your password for the `sudo` command.

```bash
user:~/washu$ cd ./docker/biolabs/uli-pipeline
user:~/washu/docker/biolabs/uli-pipeline$ sudo docker build -t biolabls/uli-pipeline .
[sudo] password for user:
```

Building the image will take some time.


Prepare data
---------

Choose a working directory.
* Make sure the directory is writable.
* Make sure there's enough disk space free (the necessary amount depends on how much data you intend to process).

Make subdirectories named `fastq` and `indexes`.

```bash
user:~/work$ mkdir fastq indexes
```

Place your FASTQ files in `fastq` subdirectory.
* If you have input files, rename them so that their name contains `input`.

```bash
user:~/work/fastq$ mv SRR1234567890.fastq SRR1234567890_input.fastq
```

If you have prebuilt bowtie indexes for your genome build, place them in `indexes/<build>`.
* If you don't have the indexes, they will be generated automatically.

```bash
user:~/work/indexes$ mkdir mm9
user:~/work/indexes$ cp ~/bowtie/mm9* ./mm9/
```

Launch pipeline
----------

Note the location of your working directory and of the cloned `washu` repository.
Run `sudo docker run -v <washu>:/washu -v <working_directory>:/data -e WASHU_ROOT=/washu -e LOCAL_USER_ID=`id -u $USER` -it biolabs/uli-pipeline /bin/bash -c "bash /washu/pipeline_chipseq.sh"`.
You might need to provide your password for the `sudo` command. The pipeline should take from several hours to several days depending on the available computing capabilities and the amount of data.

```bash
sudo docker run -v /home/user/washu:/washu -v /home/user/work:/data -e WASHU_ROOT=/washu -e LOCAL_USER_ID=`id -u $USER` -it biolabs/uli-pipeline /bin/bash -c "bash /washu/pipeline_chipseq.sh"
[sudo] password for user:
```
