import numpy as np
import pandas as pd
import pytest
from pathlib import Path
from bed.bedtrace import Bed
from scripts.util import age  # nopep8
from test.test_bed_metrics import assert_image
from pandas.util.testing import assert_frame_equal
from reports.peak_metrics import calc_consensus, venn_consensus, bar_consensus, _groups_sizes, \
    cumulative_consensus, calc_frip, frip_peaks, frip_boxplot, length_bar_plots, \
    _calculate_lengths  # nopep8
from test.fixtures import test_data, tmp_dir

threads_num = 30


@pytest.mark.parametrize("od_consensus,yd_consensus,yd_od_int", [
    ("od_consensus.bed", "yd_consensus.bed", "yd_od_int.bed")
])
def test_calc_consensus(test_data, od_consensus, yd_consensus, yd_od_int):
    expected_od_consensus = test_data("metrics/" + od_consensus)
    expected_yd_consensus = test_data("metrics/" + yd_consensus)
    expected_yd_od_int = test_data("metrics/" + yd_od_int)
    od_paths_map = {"A.bed": Bed(test_data("bed/A.bed")), "B.bed": Bed(test_data("bed/B.bed"))}
    yd_paths_map = {"C.bed": Bed(test_data("bed/C.bed")), "D.bed": Bed(test_data("bed/D.bed"))}

    od_consensus_bed, yd_consensus_bed, yd_od_int_bed = \
        calc_consensus(od_paths_map, yd_paths_map, 2.0)
    assert Bed(expected_od_consensus).cat() == od_consensus_bed.cat()
    assert Bed(expected_yd_consensus).cat() == yd_consensus_bed.cat()
    assert Bed(expected_yd_od_int).cat() == yd_od_int_bed.cat()


@pytest.mark.parametrize("fname", [
    "venn_consensus.png"
])
def test_venn_consensus(tmp_dir, test_data, fname):
    expected = test_data("metrics/" + fname)
    result = tmp_dir + "/venn_consensus.png"

    od_paths_map = {"A.bed": Bed(test_data("bed/A.bed")), "B.bed": Bed(test_data("bed/B.bed"))}
    yd_paths_map = {"C.bed": Bed(test_data("bed/C.bed")), "D.bed": Bed(test_data("bed/D.bed"))}

    od_consensus_bed, yd_consensus_bed, yd_od_int_bed = \
        calc_consensus(od_paths_map, yd_paths_map, 2.0)
    venn_consensus(od_consensus_bed, yd_consensus_bed, 2.0, result)

    assert_image(expected, result)


@pytest.mark.parametrize("fname", [
    "bar_consensus.png"
])
def test_bar_consensus(tmp_dir, test_data, fname):
    expected = test_data("metrics/" + fname)
    result = tmp_dir + "/bar_consensus.png"

    od_paths_map = {"A.bed": Bed(test_data("bed/A.bed")), "B.bed": Bed(test_data("bed/B.bed"))}
    yd_paths_map = {"C.bed": Bed(test_data("bed/C.bed")), "D.bed": Bed(test_data("bed/D.bed"))}

    od_consensus_bed, yd_consensus_bed, yd_od_int_bed = \
        calc_consensus(od_paths_map, yd_paths_map, 2.0)
    bar_consensus(od_paths_map, yd_paths_map, od_consensus_bed, yd_consensus_bed, yd_od_int_bed, 30,
                  result)

    assert_image(expected, result)


@pytest.mark.parametrize("od_consensus,yd_consensus,yd_od_int", [
    ("od_consensus.bed", "yd_consensus.bed", "yd_od_int.bed")
])
def test_groups_sizes(test_data, od_consensus, yd_consensus, yd_od_int):
    od_consensus = Bed(test_data("metrics/" + od_consensus))
    yd_consensus = Bed(test_data("metrics/" + yd_consensus))
    yd_od_int = Bed(test_data("metrics/" + yd_od_int))

    name, common, own_group, opposite_group, personal = \
        _groups_sizes(["A.bed", Bed(test_data("bed/A.bed"))], yd_od_int, od_consensus, yd_consensus)
    assert name == "A.bed"
    assert 2 == common
    assert 1 == own_group
    assert 1 == opposite_group
    assert 0 == personal


@pytest.mark.parametrize("fname", [
    "cumulative_consensus.png"
])
def test_cumulative_consensus(tmp_dir, test_data, fname):
    expected = test_data("metrics/" + fname)
    result = tmp_dir + "/cumulative_consensus.png"

    cumulative_consensus([test_data("bed/A.bed"), test_data("bed/B.bed"), test_data("bed/C.bed"),
                          test_data("bed/D.bed")], result)

    assert_image(expected, result)


@pytest.mark.parametrize("fname", [
    "df.csv"
])
def test_calc_frip(test_data, fname):
    expected_df = pd.read_csv(test_data("rip/" + fname))
    expected_df.index = [age(expected_df.loc[n]["file"]) for n in expected_df.index]
    rip_paths = [str(rip_path) for rip_path in Path(test_data("rip/")).glob("*_rip.csv")]

    age_labels, df = calc_frip(rip_paths)
    assert ['YDS', 'ODS'] == age_labels
    assert_frame_equal(expected_df, df.sort_index(), check_dtype=False)


@pytest.mark.parametrize("fname", [
    "frip_peaks.png"
])
def test_frip_peaks(tmp_dir, test_data, fname):
    expected = test_data("rip/" + fname)
    result = tmp_dir + "/frip_peaks.png"
    rip_paths = [str(rip_path) for rip_path in Path(test_data("rip/")).glob("*_rip.csv")]

    age_labels, df = calc_frip(rip_paths)
    frip_peaks(age_labels, df, result)

    assert_image(expected, result)


@pytest.mark.parametrize("fname", [
    "frip_boxplot.png"
])
def test_frip_boxplot(tmp_dir, test_data, fname):
    expected = test_data("rip/" + fname)
    result = tmp_dir + "/frip_boxplot.png"
    rip_paths = [str(rip_path) for rip_path in Path(test_data("rip/")).glob("*_rip.csv")]

    age_labels, df = calc_frip(rip_paths)
    frip_boxplot(age_labels, df, result)

    assert_image(expected, result)


@pytest.mark.parametrize("fname", [
    "length_bar_plots.png"
])
def test_length_bar_plots(tmp_dir, test_data, fname):
    expected = test_data("metrics/" + fname)
    result = tmp_dir + "/length_bar_plots.png"

    length_bar_plots([test_data("bed/A.bed"), test_data("bed/B.bed"), test_data("bed/C.bed"),
                      test_data("bed/D.bed")], 1.0, 3.0, threads_num, result)

    assert_image(expected, result)


def test_calculate_lengths(test_data):
    bins = np.logspace(1.0, 3.0, 80)

    track_path, lengths, max_bar_height = _calculate_lengths(test_data("bed/A.bed"), bins)
    assert test_data("bed/A.bed") == track_path
    assert [100] * 4 == lengths
    assert 4 == max_bar_height
