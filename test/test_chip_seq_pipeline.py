from subprocess import call
import os
import shutil
import glob

import pytest
from test.fixtures import test_data, tmp_dir


def run_pipeline():
    os.chdir(os.path.expanduser("~"))
    # os.chdir(tmp_dir)

    call(["tar", "xfz", "chip-seq-pipeline-test-data.tar.gz"])

    call(["mv", "chip-seq-pipeline-test-data", "test_data"])

    shutil.copytree("./test_data/fastq/", "./H3K4me3")

    os.makedirs("/scratch/artyomov_lab_aging/Y20O20/chipseq/indexes")

    shutil.copytree("./test_data/index/hg19",
                    "/scratch/artyomov_lab_aging/Y20O20/chipseq/indexes/hg19")

    os.mkdir("./H3K4me3_bams")

    shutil.copy("./test_data/deadzones-k36-hg19.bed",
                "./H3K4me3_bams/deadzones-k36-hg19.bed")

    call(["python", "/washu/pipeline_chipseq.py", "/root/H3K4me3"])


def check_files_not_empty(pattern, expected_files_number=None):
    files = glob.glob(pattern)

    msg = "Expected {} files for pattern '{}'".format(expected_files_number,
                                                      pattern)
    assert len(files) == expected_files_number, msg

    for f in files:
        assert os.path.getsize(f) > 0, "File {} is empty".format(f)


@pytest.mark.parametrize("ptn,expected_file_num", [
    ("./H3K4me3_bams/*.bam", 5),
    ("./H3K4me3_bams/*.bai", 5),
    ("./H3K4me3_bams_macs2_broad_0.01/*.broadPeak", 4),
    ("./H3K4me3_bams_macs2_broad_0.01/*.broadPeak_rip.csv", 4),
    ("./H3K4me3_bams_macs2_broad_0.05/*.broadPeak", 4),
    ("./H3K4me3_bams_macs2_broad_0.05/*.broadPeak_rip.csv", 4),
    ("./H3K4me3_bams_macs2_broad_0.1/*.broadPeak", 4),
    ("./H3K4me3_bams_macs2_broad_0.1/*.broadPeak_rip.csv", 4),
    ("./H3K4me3_bams_macs2_q0.01/*.narrowPeak", 4),
    ("./H3K4me3_bams_macs2_q0.01/*.narrowPeak_rip.csv", 4),
    ("./H3K4me3_bams_macs2_q0.05/*.narrowPeak", 4),
    ("./H3K4me3_bams_macs2_q0.05/*.narrowPeak_rip.csv", 4),
    ("./H3K4me3_bams_macs2_q0.1/*.narrowPeak", 4),
    ("./H3K4me3_bams_macs2_q0.1/*.narrowPeak_rip.csv", 4),
    ("./H3K4me3_bams_rpkms/*.bw", 5),
    ("./H3K4me3_bams_rseg/*_domains.bed", 4),
    ("./H3K4me3_bams_sicer_0.01/*-island.bed", 4),
    ("./H3K4me3_bams_sicer_0.01/*-island.bed_rip.csv", 4),
    ("./H3K4me3_bams_unique/*.bam", 5)
])
def check_files(ptn, expected_file_num):
    os.chdir(os.path.expanduser("~"))
    check_files_not_empty(ptn, expected_file_num)


def test_pipeline():
    run_pipeline()
    check_files()
