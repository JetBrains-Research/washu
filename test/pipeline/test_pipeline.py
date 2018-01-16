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
                            line, flags=re.IGNORECASE) is not None and \
                                re.match('.*ValueError.*', line) is None:
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
    check_files("fastq_bams_unique_tags_bws/*.bw", 6)


def test_rpkm_step():
    check_files("fastq_bams_rpkms/*.bw", 6)


def test_macs2_broad():
    check_files("fastq_bams_macs2_broad_0.1/*.broadPeak", 4)
    check_files("fastq_bams_macs2_broad_0.1/*.broadPeak_rip.csv", 4)
    check_files("fastq_bams_macs2_broad_0.1/macs2_report.csv", 1)


def test_macs2_narrow():
    check_files("fastq_bams_macs2_q0.01/*.narrowPeak", 4)
    check_files("fastq_bams_macs2_q0.01/*.narrowPeak_rip.csv", 4)
    check_files("fastq_bams_macs2_q0.01/macs2_report.csv", 1)


def test_rseg():
    check_files("fastq_bams_rseg/*_domains.bed", 4)
    check_files("fastq_bams_rseg/peaks_report.csv", 1)


def test_sicer():
    check_files("fastq_bams_sicer/*island.bed", 4)
    check_files("fastq_bams_sicer/*island.bed_rip.csv", 4)
    check_files("fastq_bams_sicer/*removed-1.bed", 0)
    # BATCH is true, these shouldn't be cleaned up
    check_files("fastq_bams_sicer/*input*.bed", 4)
    check_files("fastq_bams_sicer/*pileup.bed", 6)
    check_files("fastq_bams_sicer/peaks_report.csv", 1)


def test_bam2tags():
    with tempfile.NamedTemporaryFile(
            mode='w', suffix='.tags', prefix='bam', delete=False) as tmpfile:
        run([["bash", "/washu/scripts/bam2tags.sh",
              os.path.expanduser("~/fastq_bams/OD1_k4me3_hg19.bam"), "150"],
             ["head"]], stdout=tmpfile)
        assert Path(tmpfile.name).read_text() == Path("data/bam2tags.tag").read_text()


def test_reads2bam():
    bam = run([["bash", "/washu/scripts/reads2bam.sh",
                os.path.expanduser("~/data/reads.bed"),
                os.path.expanduser("~/index/hg19/hg19.chrom.sizes")]])[0]. \
        decode('utf-8').strip()
    assert os.path.expanduser("~/data/reads.bam") == bam
    reads1 = Path("data/reads.bed").read_text()
    reads2 = run([["bedtools", "bamtobed", "-i", bam]])[0].decode('utf-8')
    assert reads1 == reads2


def test_tags2bdg():
    bdg = run([["bash", "/washu/scripts/tags2bdg.sh",
                os.path.expanduser("~/data/bam2tags.tag")]])[0].decode('utf-8')
    assert Path("data/bam2tags.bdg").read_text() == bdg


def test_signals():
    if not os.path.exists("regions.bed"):
        shutil.copy("washu_test_data/data/regions.bed", "regions.bed")
    if os.path.exists("fastq_bams_bws/regions"):
        shutil.rmtree("fastq_bams_bws/regions")
    call(["bash", "/washu/downstream/signals/signals_bw.sh",
          os.path.expanduser("~/fastq_bams_unique_tags_bws"),
          os.path.expanduser("~/data/regions.bed"),
          os.path.expanduser("~/index/hg19/hg19.chrom.sizes")])

    check_files("fastq_bams_unique_tags_bws/hg19.chrom.sizes.tsv", 1)
    check_files("fastq_bams_unique_tags_bws/regions/*_signal.log", 1)
    check_files("fastq_bams_unique_tags_bws/regions/regions.tsv", 1)
    check_files("fastq_bams_unique_tags_bws/regions/regions_raw.tsv", 1)
    check_files("fastq_bams_unique_tags_bws/regions/regions_raw_data.csv", 1)
    check_files("fastq_bams_unique_tags_bws/regions/regions_raw.png", 1)
    check_files("fastq_bams_unique_tags_bws/regions/regions_raw_pca.png", 1)
    check_files("fastq_bams_unique_tags_bws/regions/regions_raw_pca_fit_error.csv", 1)
    check_files("fastq_bams_unique_tags_bws/regions/regions_raw_t.tsv", 1)
    check_files("fastq_bams_unique_tags_bws/regions/regions_raw_u.tsv", 1)
    check_files("fastq_bams_unique_tags_bws/regions/regions_q.tsv", 1)
    check_files("fastq_bams_unique_tags_bws/regions/regions_rpm.tsv", 1)
    check_files("fastq_bams_unique_tags_bws/regions/regions_rpkm.tsv", 1)
    check_files("fastq_bams_unique_tags_bws/regions/regions_scores.tsv", 1)
    # TODO[shpynov] fix
    # check_files("fastq_bams_unique_tags_bws/regions/regions_tmm_scores.tsv", 1)
