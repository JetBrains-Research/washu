import tempfile
from subprocess import call
import os
import shutil
import glob

import pytest

from pathlib import Path
from scripts.util import run
from test.fixtures import test_data, tmp_dir


def setup_module(module):
    """ setup any state specific to the execution of the given module."""
    os.chdir(os.path.expanduser("~"))
    if not os.path.exists("./washu_test_data"):
        call(["tar", "xfz", "washu_test_data.tar.gz"])
        if os.path.exists("./fastq"):
            shutil.rmtree("./fastq")
        shutil.copytree("./washu_test_data/fastq/", "./fastq")
        if os.path.exists("./index"):
            shutil.rmtree("./index")
        shutil.copytree("./washu_test_data/index", "./index")
        if os.path.exists("./data"):
            shutil.rmtree("./data")
        shutil.copytree("./washu_test_data/data", "./data")

    run_pipeline()


def run_pipeline():
    # Unpack test data and prepare layout
    os.chdir(os.path.expanduser("~"))
    if os.path.exists("./pipeline_finished.txt"):
        return

    # Prepare regions for RSEG
    if not os.path.exists("./fastq_bams"):
        os.mkdir("./fastq_bams")
        if not os.path.exists("./fastq_bams/deadzones-k36-hg19.bed"):
            shutil.copy("./index/hg19/deadzones-k36-hg19.bed",
                        "./fastq_bams/deadzones-k36-hg19.bed")

    # Launch pipeline
    call(["python", "/washu/pipeline_chipseq.py", "/root/fastq", "/root/index", "hg19"])

    # Touch marker file
    open("./pipeline_finished.txt", "a")


def check_files(pattern, expected_files_number=None):
    files = glob.glob(pattern)
    msg = "Expected {} files for pattern '{}'".format(expected_files_number,
                                                      pattern)
    assert len(files) == expected_files_number, msg
    for f in files:
        assert os.path.getsize(f) > 0, "File {} is empty".format(f)


def test_pipeline_bams():
    '''NOTE: Pipeline tests are split into different tests for better problems reporting!'''
    check_files("./fastq_bams/*.bam", 6)


def test_pipeline_bw():
    check_files("./fastq_bams_bws/*.bw", 6)


def test_pipeline_unique():
    check_files("./fastq_bams_unique/*.bam", 6)


def test_pipeline_unique_tags_bws():
    check_files("./fastq_bams_unique_tags_bws/*.bw", 6)


def test_pipeline_rpkm_step():
    check_files("./fastq_bams_rpkms/*.bw", 6)


def test_pipeline_macs2_broad():
    # Peak calling steps with different settings MACS2 --broad
    check_files("./fastq_bams_macs2_broad_0.1/*.broadPeak", 4)
    check_files("./fastq_bams_macs2_broad_0.1/*.broadPeak_rip.csv", 4)
    check_files("./fastq_bams_macs2_broad_0.05/*.broadPeak", 4)
    check_files("./fastq_bams_macs2_broad_0.05/*.broadPeak_rip.csv", 4)
    check_files("./fastq_bams_macs2_broad_0.01/*.broadPeak", 4)
    check_files("./fastq_bams_macs2_broad_0.01/*.broadPeak_rip.csv", 4)


def test_pipeline_macs2_narrow():
    # Peak calling steps with different settings MACS2 narrow
    check_files("./fastq_bams_macs2_q0.1/*.narrowPeak", 4)
    check_files("./fastq_bams_macs2_q0.1/*.narrowPeak_rip.csv", 4)
    check_files("./fastq_bams_macs2_q0.05/*.narrowPeak", 4)
    check_files("./fastq_bams_macs2_q0.05/*.narrowPeak_rip.csv", 4)
    check_files("./fastq_bams_macs2_q0.01/*.narrowPeak", 4)
    check_files("./fastq_bams_macs2_q0.01/*.narrowPeak_rip.csv", 4)


def test_pipeline_rseg():
    check_files("./fastq_bams_rseg/*_domains.bed", 4)


def test_pipeline_sicer():
    check_files("./fastq_bams_sicer_0.01/*-island.bed", 4)
    check_files("./fastq_bams_sicer_0.01/*-island.bed_rip.csv", 4)


def test_bam2tags():
    with tempfile.NamedTemporaryFile(
            mode='w', suffix='.tags', prefix='bam', delete=False) as tmpfile:
        run([["bash", "/washu/scripts/bam2tags.sh", "/root/fastq_bams/OD1_k4me3_hg19.bam", "150"],
             ["head"]], stdout=tmpfile)
        assert Path(tmpfile.name).read_text() == Path("/root/data/bam2tags.tag").read_text()


def test_reads2bam():
    bam = run([["bash", "/washu/scripts/reads2bam.sh",
                "/root/data/reads.bed", "/root/index/hg19/hg19.chrom.sizes"]])[0]. \
        decode('utf-8').strip()
    assert "/root/data/reads.bam" == bam
    reads1 = Path("/root/data/reads.bed").read_text()
    reads2 = run([["bedtools", "bamtobed", "-i", bam]])[0].decode('utf-8')
    assert reads1 == reads2


def test_tags2bdg():
    bdg = run([["bash", "/washu/scripts/tags2bdg.sh", "/root/data/bam2tags.tag"]])[0]. \
        decode('utf-8')
    assert Path("./data/bam2tags.bdg").read_text() == bdg


def test_signals():
    if not os.path.exists("./regions.bed"):
        shutil.copy("./washu_test_data/data/regions.bed", "./regions.bed")
    if os.path.exists("./fastq_bams_bws/regions"):
        shutil.rmtree("./fastq_bams_bws/regions")
    call(["bash", "/washu/parallel/signals_bw.sh",
          "/root/fastq_bams_unique_tags_bws", "/root/data/regions.bed", "regions",
          "/root/index/hg19/hg19.chrom.sizes"])

    os.chdir(os.path.expanduser("~"))
    check_files("./fastq_bams_unique_tags_bws/hg19.chrom.sizes.tsv", 1)
    check_files("./fastq_bams_unique_tags_bws/regions/*_signal.log", 1)
    check_files("./fastq_bams_unique_tags_bws/regions/regions.tsv", 1)
    check_files("./fastq_bams_unique_tags_bws/regions/regions_raw.tsv", 1)
    check_files("./fastq_bams_unique_tags_bws/regions/regions.tsv", 1)
    assert Path("./fastq_bams_unique_tags_bws/regions/regions_raw.tsv").read_text() == \
           Path("./data/regions_raw.tsv").read_text()  # nopep8
