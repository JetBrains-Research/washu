Docker Image with test data
===========================

This is just an image `continuumio/miniconda` with test data from
`/mnt/stripe/washu_test_data`.

Build
-----
To build you have to copy file `washu_test_data.tar.gz` to folder with this image.

```bash
# Copy test data to folder with Docker file
tar -cvzf washu_test_data.tar.gz -C /mnt/stripe/ washu_test_data

# Build Docker
docker build -t biolabs/test-data .

# Clean up
rm washu_test_data.tar.gz
```

Push
----
Before push you have to login to docker hub first.
```bash
docker login -u biolabs
```

Then you just push current image 
```bash
docker push biolabs/test-data
```

Test data
---------

* `data` - test data for random scripts in tools in `/scripts` folder
* `fastq` - 4 k4me3 fastq files limited to chr22
* `index` - indices for hg19 genome limiter to chr22
```
> tree washu_test_data
washu_test_data/
├── data
│   ├── bam2tags.bdg
│   ├── bam2tags.tag
│   ├── reads.bed
│   ├── regions.bed
│   └── regions_raw.tsv
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