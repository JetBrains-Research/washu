from subprocess import call
import os
import shutil
import glob

import pytest
from test.fixtures import test_data, tmp_dir


def check_files(pattern, expected_files_number=None):
    files = glob.glob(pattern)
    msg = "Expected {} files for pattern '{}'".format(expected_files_number,
                                                      pattern)
    assert len(files) == expected_files_number, msg
    for f in files:
        assert os.path.getsize(f) > 0, "File {} is empty".format(f)


def run_pipeline():
    # Unpack test data and prepare layout
    os.chdir(os.path.expanduser("~"))
    call(["tar", "xfz", "chip-seq-pipeline-test-data.tar.gz"])
    call(["mv", "chip-seq-pipeline-test-data", "test_data"])
    shutil.copytree("./test_data/fastq/", "./H3K4me3")
    os.makedirs("/scratch/artyomov_lab_aging/Y20O20/chipseq/indexes")
    shutil.copytree("./test_data/index/hg19",
                    "/scratch/artyomov_lab_aging/Y20O20/chipseq/indexes/hg19")
    os.mkdir("./H3K4me3_bams")
    shutil.copy("./test_data/deadzones-k36-hg19.bed",
                "./H3K4me3_bams/deadzones-k36-hg19.bed")
    # Launch pipeline
    call(["python", "/washu/pipeline_chipseq.py", "/root/H3K4me3"])


def check_all():
    os.chdir(os.path.expanduser("~"))
    check_files("./H3K4me3_bams/*.bam", 5)
    check_files("./H3K4me3_bams/*.bai", 5)

    check_files("./H3K4me3_bams_macs2_broad_0.01/*.broadPeak", 4)
    check_files("./H3K4me3_bams_macs2_broad_0.01/*.broadPeak_rip.csv", 4)

    check_files("./H3K4me3_bams_macs2_broad_0.05/*.broadPeak", 4)
    check_files("./H3K4me3_bams_macs2_broad_0.05/*.broadPeak_rip.csv", 4)

    check_files("./H3K4me3_bams_macs2_broad_0.1/*.broadPeak", 4)
    check_files("./H3K4me3_bams_macs2_broad_0.1/*.broadPeak_rip.csv", 4)

    check_files("./H3K4me3_bams_macs2_q0.01/*.narrowPeak", 4)
    check_files("./H3K4me3_bams_macs2_q0.01/*.narrowPeak_rip.csv", 4)

    check_files("./H3K4me3_bams_macs2_q0.05/*.narrowPeak", 4)
    check_files("./H3K4me3_bams_macs2_q0.05/*.narrowPeak_rip.csv", 4)

    check_files("./H3K4me3_bams_macs2_q0.1/*.narrowPeak", 4)
    check_files("./H3K4me3_bams_macs2_q0.1/*.narrowPeak_rip.csv", 4)

    check_files("./H3K4me3_bams_rpkms/*.bw", 5)

    check_files("./H3K4me3_bams_rseg/*_domains.bed", 4)

    check_files("./H3K4me3_bams_sicer_0.01/*-island.bed", 4)
    check_files("./H3K4me3_bams_sicer_0.01/*-island.bed_rip.csv", 4)

    check_files("./H3K4me3_bams_unique/*.bam", 5)


def run_signals_bw():
    os.chdir(os.path.expanduser("~"))
    shutil.copytree("./test_data/regions.bed", "./regions.bed")
    call(["bash",  "/washu/parallel/signals_bw.py",
          "/root/H3K4me3_bams_bws", "/root/regions.bed", "regions",
         "/scratch/artyomov_lab_aging/Y20O20/chipseq/indexes/hg19/hg19.chrom.sizes"])


def check_signals_bw():
    os.chdir(os.path.expanduser("~"))
    check_files("./H3K4me3_bams_bws/regions/*_pca.png", 5)
    check_files("./H3K4me3_bams_bws/regions/*_fit_error.csv", 5)


def test_pipeline():
    run_pipeline()
    run_signals_bw()
    check_all()
    check_signals_bw()
