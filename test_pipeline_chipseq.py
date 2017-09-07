#!/usr/bin/env python

"""
Test for a chip-seq technical pipeline.

To run test run Docker with mounted test_data folder.

docker run -v /mnt/stripe/washu:/washu -v /mnt/stripe/chip-seq-pipeline-test-data:/root/test_data -t -i washu /bin/bash

In Docker run command.

source activate py3.5 && python /washu/test_pipeline_chipseq.py

"""

from subprocess import call
import os
import shutil
import glob


def run_pipeline():
    os.chdir(os.path.expanduser("~"))

    call(["tar", "xfz", "chip-seq-pipeline-test-data.tar.gz"])

    call(["mv", "chip-seq-pipeline-test-data", "test_data"])

    shutil.copytree("./test_data/fastq/", "./H3K4me3")

    os.makedirs("/scratch/artyomov_lab_aging/Y20O20/chipseq/indexes")

    shutil.copytree("./test_data/index/hg19", "/scratch/artyomov_lab_aging/Y20O20/chipseq/indexes/hg19")

    os.mkdir("./H3K4me3_bams")

    shutil.copy("./test_data/deadzones-k36-hg19.bed", "./H3K4me3_bams/deadzones-k36-hg19.bed")

    call(["python", "/washu/pipeline_chipseq.py", "/root/H3K4me3"])


def test_files_not_empty(pattern, expected_files_number=None):
    files = glob.glob(pattern)
    if expected_files_number and len(files) != expected_files_number:
        print("Expected {} files for pattern '{}', but {} found.".format(expected_files_number, pattern, len(files)))
        exit(1)
    for f in files:
        if os.path.getsize(f) == 0:
            print("File {} is empty.".format(f))
            exit(1)


def test_files():
    os.chdir(os.path.expanduser("~"))
    test_files_not_empty("./H3K4me3_bams/*.bam", 5)
    test_files_not_empty("./H3K4me3_bams/*.bai", 5)

    test_files_not_empty("./H3K4me3_bams_macs2_broad_0.01/*.broadPeak", 4)
    test_files_not_empty("./H3K4me3_bams_macs2_broad_0.01/*.broadPeak_rip.csv", 4)

    test_files_not_empty("./H3K4me3_bams_macs2_broad_0.05/*.broadPeak", 4)
    test_files_not_empty("./H3K4me3_bams_macs2_broad_0.05/*.broadPeak_rip.csv", 4)

    test_files_not_empty("./H3K4me3_bams_macs2_broad_0.1/*.broadPeak", 4)
    test_files_not_empty("./H3K4me3_bams_macs2_broad_0.1/*.broadPeak_rip.csv", 4)

    test_files_not_empty("./H3K4me3_bams_macs2_q0.01/*.narrowPeak", 4)
    test_files_not_empty("./H3K4me3_bams_macs2_q0.01/*.narrowPeak_rip.csv", 4)

    test_files_not_empty("./H3K4me3_bams_macs2_q0.05/*.narrowPeak", 4)
    test_files_not_empty("./H3K4me3_bams_macs2_q0.05/*.narrowPeak_rip.csv", 4)

    test_files_not_empty("./H3K4me3_bams_macs2_q0.1/*.narrowPeak", 4)
    test_files_not_empty("./H3K4me3_bams_macs2_q0.1/*.narrowPeak_rip.csv", 4)

    test_files_not_empty("./H3K4me3_bams_rpkms/*.bw", 5)

    test_files_not_empty("./H3K4me3_bams_rseg/*_domains.bed", 4)

    test_files_not_empty("./H3K4me3_bams_sicer_0.01/*-island.bed", 4)
    test_files_not_empty("./H3K4me3_bams_sicer_0.01/*-island.bed_rip.csv", 4)

    test_files_not_empty("./H3K4me3_bams_unique/*.bam", 5)

if __name__ == '__main__':
    run_pipeline()
    test_files()
