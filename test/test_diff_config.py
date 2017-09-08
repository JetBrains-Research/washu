from bed.bedtrace import run
from pipeline_utils import run_bash, PROJECT_ROOT_PATH

import pytest
from test.fixtures import test_data


def test_diffbind_config(capfd, test_data):
    run_bash("scripts/diff_config.sh",
             test_data("diff_config/k4me3_bams"),
             test_data("diff_config/k4me3_bams_macs2_broad_0.1"))

    out, _err = capfd.readouterr()
    res = out.replace(test_data("diff_config/"), "").replace(PROJECT_ROOT_PATH,
                                                             ".")
    print(res)
    assert res == """bash ./scripts/diff_config.sh k4me3_bams \
k4me3_bams_macs2_broad_0.1
SampleID,Tissue,Factor,Condition,Replicate,bamReads,ControlID,bamControl,Peaks,\
PeakCaller
OD1,CD14,Age,O,1,k4me3_bams/OD1_k4me3_hg19.bam,O_input,OD_input.bam,\
k4me3_bams_macs2_broad_0.1/OD1_k4me3_hg19_broad_0.1_peaks.xls,macs
"""
