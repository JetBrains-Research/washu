Docker Image with test data and tools
=====================================

This is just an image `ubuntu:latest` - Ubuntu LTS with 
* all the environment setup
* test data from `/mnt/stripe/washu_test_data` 

Prepare Test data
---------
```bash
# Create .tar.gz with test data
tar -cvzf washu_test_data.tar.gz -C /mnt/stripe/washu_test_data
```
* `fastq` - fastq files (chr22)
* `index` - indices for hg19 genome (chr22)
```
> tree washu_test_data
washu_test_data/
├── fastq
│   ├── OD1_k4me3.fq
│   ├── OD3_k4me3.fq
│   ├── od_input.fq
│   ├── YD1_k4me3.fq
│   ├── YD3_k4me3.fq
│   └── yd_input.fq
└── index
    └── hg19
        ├── chr22.fa
        ├── deadzones-k36-hg19.bed
        └── hg19.chrom.sizes
```

Build
-----
Change directory to `$WASHU_ROOT/docker/biolabs/washu` and copy `washu_test_data.tar.gz` file
 to current directory.

```bash
docker build -t biolabs/washu .
```

Push
----
Before push you have to login to docker hub first.
```bash
docker login -u biolabs
```

Then you just push current image 
```bash
docker push biolabs/washu
```