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


def test_bw():
    check_files("fastq_bams_bws/*.bw", 6)


def test_unique():
    check_files("fastq_bams_unique/*.bam", 6)


def test_bam_qc():
    check_files("fastq_bams/qc/*.pbc_nrf.tsv", 6)
    # TODO[shpynov] fix bam_qc.sh: could not find function "runmean"
    # check_files("fastq_bams/qc/*.phantom.tsv", 6)
    # check_files("fastq_bams/qc/*.pdf", 6)


def test_macs2_broad():
    check_files("fastq_bams_macs2_broad_0.1/*.broadPeak", 4)
    check_files("fastq_bams/pileup/*pileup.bed", 6)


def test_macs2_narrow():
    check_files("fastq_bams_macs2_q0.05/*.narrowPeak", 4)
    check_files("fastq_bams/pileup/*pileup.bed", 6)


def test_rseg():
    check_files("fastq_bams_rseg/*_domains.bed", 4)
    check_files("fastq_bams/pileup/*pileup.bed", 6)


def test_sicer():
    check_files("fastq_bams/pileup/*pileup.bed", 6)
    check_files("fastq_bams_sicer/*island.bed", 4)
    check_files("fastq_bams_sicer/*removed-1.bed", 0)
    check_files("fastq_bams_sicer/*input*.bed", 2)


def test_span():
    check_files("fastq_bams_span/fit/*.span", 4)
    check_files("fastq_bams_span/logs/*.log", 4)
    check_files("fastq_bams_span/*_peaks.bed", 4)


def test_reads2bam():
    bam = run([["bash", "/washu/scripts/reads2bam.sh",
                "/washu/test/testdata/pileup/a_pileup.bed",
                os.path.expanduser("~/index/hg19/hg19.chrom.sizes")]])[0]. \
        decode('utf-8').strip()
    assert "/washu/test/testdata/pileup/a_pileup.bam" == bam
    reads1 = Path("/washu/test/testdata/pileup/a_pileup.bed").read_text()
    reads2 = run([["bedtools", "bamtobed", "-i", bam]])[0].decode('utf-8')
    # Ignore trailing 0
    assert re.sub('\\t0', '', reads1).strip() == reads2.strip()


def test_signals():
    # Create signals folder
    signals_path = os.path.expanduser("~/signals/H3K4me3")
    if not os.path.exists(signals_path):
        os.makedirs(signals_path)

    # Create symbolic links for BAMs
    bams_path = os.path.expanduser("~/fastq_bams")
    for r in glob.glob("{}/*.bam".format(bams_path)):
        call(["bash", "-c", "ln -sf {} {}/".format(r, signals_path)])

    # Process single file
    call(["bash", "/washu/downstream/signals/signals.sh",
          signals_path,
          "120",
          os.path.expanduser("/washu/test/testdata/signal/regions.bed"),
          os.path.expanduser("~/index/hg19/hg19.chrom.sizes"),
          os.path.expanduser("/washu/test/testdata/bed/peaks.bed")])

    # And check
    check_files(signals_path + "/120/*.bw", 6)
    check_files(signals_path + "/120/tags_bw_logs/*.log", 6)
    check_files(signals_path + "/120/hg19.chrom.sizes.tsv", 1)
    check_files(signals_path + "/120/peaks.bed.tsv", 1)
    check_files(signals_path + "/120/*_tsv.log", 1)
    check_files(signals_path + "/120/regions/*_signal.log", 1)
    check_files(signals_path + "/120/regions/regions.tsv", 1)
    check_files(signals_path + "/120/regions/regions_raw.tsv", 1)
    check_files(signals_path + "/120/regions/*.tsv", 8)
    check_files(signals_path + "/120/regions/*.png", 7)
    check_files(signals_path + "/120/regions/*_pca_fit_error.csv", 7)

    # Process folder
    call(["bash", "/washu/downstream/signals/signals.sh",
          signals_path,
          "150",
          os.path.expanduser("/washu/test/testdata/signal"),
          os.path.expanduser("~/index/hg19/hg19.chrom.sizes"),
          os.path.expanduser("/washu/test/testdata/bed/peaks.bed")])

    # And check
    check_files(signals_path + "/150/a/a1/a1*.tsv", 8)
    check_files(signals_path + "/150/a/a2/a2*.tsv", 8)
    check_files(signals_path + "/150/b/c/c/c*.tsv", 8)
    check_files(signals_path + "/150/regions/regions*.tsv", 8)

    # Create and check report
    report_path = os.path.expanduser("~/signals/report.tsv")
    with open(report_path, 'w') as r:
        call(["bash", "/washu/downstream/signals/signals_report.sh",
              os.path.expanduser("~/signals")], stdout=r)
    check_files(report_path, 1)
    lines = {line.rstrip('\n') for line in open(report_path)}
    assert 'H3K4me3\t150\tb/c/c\tfripm_1kbp\t0' in lines
    assert 'H3K4me3\t150\ta/a1\trawq\t0' in lines
    assert 'H3K4me3\t120\tregions\traw\t1' in lines
    assert len(lines) == 36
