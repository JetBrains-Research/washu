import re
import shutil
import tempfile
from subprocess import call
import os
import glob

import pytest

from pathlib import Path
from scripts.util import run
from test.fixtures import test_data, tmp_dir


def setup_module(module):
    """ setup any state specific to the execution of the given module."""
    os.chdir(os.path.expanduser("~"))
    if os.path.exists("pipeline_finished.txt"):
        return

    # Launch pipeline
    call(["python", "/washu/pipeline_chipseq.py",
          os.path.expanduser("~/fastq"), os.path.expanduser("~/index"), "hg19"])

    # Touch marker file
    open("pipeline_finished.txt", "a")


def check_files(pattern, expected_files_number=None):
    files = glob.glob(pattern)
    msg = "Expected {} files for pattern '{}'".format(expected_files_number,
                                                      pattern)
    assert len(files) == expected_files_number, msg
    for f in files:
        assert os.path.getsize(f) > 0, "File {} is empty".format(f)


def test_errors():
    for file in glob.glob('fastq**/*.log', recursive=True):
        with open(file, 'r') as log:
            for line in log:
                # IGNORE macs2 ValueError
                # See for details: https://github.com/JetBrains-Research/washu/issues/14
                if re.match('.*(error|exception|No such file or directory).*',
                            line, flags=re.IGNORECASE) and \
                        not re.match('.*ValueError.*', line):
                    print('Error in log file', file, ':', line)
                    assert False


def test_bams():
    """NOTE: Pipeline tests are split into different tests for better problems reporting!"""
    check_files("fastq_bams/*.bam", 6)
    check_files("fastq_bams/bowtie_report.csv", 1)


def test_bw():
    check_files("fastq_bams_bws/*.bw", 6)


def test_unique():
    check_files("fastq_bams_unique/*.bam", 6)


def test_bam_qc():
    check_files("fastq_bams/qc/*.pbc_nrf.tsv", 6)
    check_files("fastq_bams/qc/*.phantom.tsv", 6)
    check_files("fastq_bams/qc/*.pdf", 6)


def test_unique_tags_bws():
    check_files("signals/*.bw", 6)


def test_rpkm_step():
    check_files("fastq_bams_rpkms/*.bw", 6)


def test_macs2_broad():
    check_files("fastq_bams_macs2_broad_0.1/*.broadPeak", 4)
    check_files("fastq_bams_macs2_broad_0.1/*.broadPeak_rip.csv", 4)
    check_files("fastq_bams_macs2_broad_0.1/macs2_report.csv", 1)
    check_files("fastq_bams/pileup/*pileup.bed", 6)


def test_macs2_narrow():
    check_files("fastq_bams_macs2_q0.05/*.narrowPeak", 4)
    check_files("fastq_bams_macs2_q0.05/*.narrowPeak_rip.csv", 4)
    check_files("fastq_bams_macs2_q0.05/macs2_report.csv", 1)
    check_files("fastq_bams/pileup/*pileup.bed", 6)


def test_rseg():
    check_files("fastq_bams_rseg/*_domains.bed", 4)
    check_files("fastq_bams_rseg/peaks_report.csv", 1)
    check_files("fastq_bams/pileup/*pileup.bed", 6)


def test_sicer():
    check_files("fastq_bams/pileup/*pileup.bed", 6)
    check_files("fastq_bams_sicer/*island.bed", 4)
    check_files("fastq_bams_sicer/*island.bed_rip.csv", 4)
    check_files("fastq_bams_sicer/*removed-1.bed", 0)
    check_files("fastq_bams_sicer/*input*.bed", 2)
    check_files("fastq_bams_sicer/peaks_report.csv", 1)


def test_reads2bam():
    bam = run([["bash", "/washu/scripts/reads2bam.sh",
                os.path.expanduser("~/data/reads.bed"),
                os.path.expanduser("~/index/hg19/hg19.chrom.sizes")]])[0]. \
        decode('utf-8').strip()
    assert os.path.expanduser("~/data/reads.bam") == bam
    reads1 = Path("data/reads.bed").read_text()
    reads2 = run([["bedtools", "bamtobed", "-i", bam]])[0].decode('utf-8')
    assert reads1 == reads2


def test_signals():
    # Copy regions.bed
    if not os.path.exists("regions.bed"):
        shutil.copy("washu_test_data/data/regions.bed", "regions.bed")
    if os.path.exists("fastq_bams_bws/regions"):
        shutil.rmtree("fastq_bams_bws/regions")

    # Create signals folder
    signals_path = os.path.expanduser("~/signals")
    if not os.path.exists(signals_path):
        os.mkdir(signals_path)

    # Create symbolic links for BAMs
    bams_path = os.path.expanduser("~/fastq_bams")
    for f in glob.glob("{}/*.bam".format(bams_path)):
        call(["bash", "-c", "ln -sf {} {}/".format(f, signals_path)])

    # Call
    call(["bash", "/washu/downstream/signals/signals.sh",
          signals_path,
          "150",
          os.path.expanduser("~/data/regions.bed"),
          os.path.expanduser("~/index/hg19/hg19.chrom.sizes"),
          os.path.expanduser("~/data/regions.bed")])

    # And check
    check_files("signals/150/hg19.chrom.sizes.tsv", 1)
    check_files("signals/150/regions/*_signal.log", 1)
    check_files("signals/150/regions/bw_signals_*.log", 1)
    check_files("signals/150/regions/regions.tsv", 1)
    check_files("signals/150/regions/regions_raw.tsv", 1)
