import os
import subprocess

import pytest
from test.fixtures import test_data

from scripts.util import find_input
from scripts.util import lcs
from scripts.util import macs_species
from pipeline_utils import PROJECT_ROOT_PATH


def test_lcs():
    lcs1 = lcs('UW_CD14_RO01746_k4me3_1_1_ENCFF001FYS.fastq',
               'Broad_CD14_2_input_ENCFF000CCW.fastq')
    lcs2 = lcs('UW_CD14_RO01746_k4me3_1_1_ENCFF001FYS.fastq',
               'UW_CD14_input_ENCFF001HUV.fastq')
    assert lcs1 < lcs2


@pytest.mark.parametrize("input,donor", [
    ("", "40_donor6_input.bam"),
    ("40_donor6_input.bam", "37_donor6_k27ac.bam"),
    ("44_donor7_input.bam", "41_donor7_k27ac.bam"),
    ("jcl320_wt1_gm_input_ctrl.1919_8.R1_mm10.bam",
     "jcl320_ko_gm_h3k27ac_chipd_dna.1919_8.R1_mm10.bam"),
    ("jcl320_wt1_gm_input_ctrl.1919_8.R1_mm10.bam",
     "jcl320_wt1_gm_bhlhe40_chipd_dna.1919_8.R1_mm10.bam"),
])
def test_find_input(test_data, input, donor):
    assert find_input(test_data("input/" + donor)) == input


@pytest.mark.parametrize("build,species", [
    ("hg18", "hs"),
    ("hg19", "hs"),
    ("hg37", "hs"),
    ("mm9", "mm"),
    ("mm10", "mm"),
])
def test_macs_input(build, species):
    assert macs_species(build) == species


@pytest.mark.parametrize("path,expected", [
    ("geo", "geo"),
    ("geo/.", "geo"),
    ("geo/tmp/..", "geo"),
    ("geo/../geo/tmp/doo/..", "geo/tmp"),
    ("geo/../geo/tmp/doo/file.txt", "geo/tmp/doo/file.txt"),
    ("geo/symlink", "geo/tmp"),
    ("geo/symlink/..", "geo"),
    ("geo/../geo/symlink/doo/file.txt", "geo/tmp/doo/file.txt"),
])
def test_expand_path(tmpdir, capfd, path, expected):
    os.makedirs(os.path.join(tmpdir, "geo/tmp/doo"))
    os.symlink(os.path.join(tmpdir, "geo/tmp"),
               os.path.join(tmpdir, "geo/symlink"))
    open(os.path.join(tmpdir, "geo/tmp/doo/file.txt"), 'a').close()

    path = os.path.join(tmpdir, path)
    util_sh = os.path.join(PROJECT_ROOT_PATH, "parallel/util.sh")
    subprocess.run(
        "bash -c 'source {}; echo $(expand_path \"{}\")'".format(util_sh,
                                                                 path),
        shell=True, check=True)
    out, _err = capfd.readouterr()
    assert out == os.path.join(tmpdir, expected) + "\n"
