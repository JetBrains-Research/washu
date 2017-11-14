import platform
import os
import sys
import shutil
from pathlib import Path
import pandas as pd
import pytest
from test.fixtures import test_data, tmp_dir


from reports.bed_metrics import bed_metric_table, color_annotator_age, \
    plot_metric_heatmap, _run_metric_jaccard, _run_metric_intersection, \
    label_converter_donor_and_tool, color_annotator_chain, \
    color_annotator_outlier  # nopep8


@pytest.mark.parametrize("jaccard,fname,swap_ab", [
    (False, "metric1.csv", False),
    (False, "metric1_ba.csv", True),
    (True, "metric_j.csv", False),
    (True, "metric_j_ba.csv", True)
])
def test_bed_metric_table(test_data, jaccard, fname, swap_ab):
    a_paths \
        = [Path(test_data("metrics/c{}.bed".format(i))) for i in range(1, 4)]
    b_paths = [Path(test_data("metrics/a.bed")), *a_paths]

    df = bed_metric_table(a_paths if not swap_ab else b_paths,
                          b_paths if not swap_ab else a_paths,
                          jaccard=jaccard)
    expected = pd.DataFrame.from_csv(test_data("metrics/{}".format(fname)))

    assert str(expected) == str(df)


@pytest.mark.parametrize("jaccard,a,b,value", [
    (False, "zero.bed", "a.bed", 0),
    (False, "a.bed", "zero.bed", 0),
    (False, "zero.bed", "zero.bed", 0),
    (True, "zero.bed", "a.bed", 0.0),
    (True, "a.bed", "zero.bed", 0.0),
    (True, "zero.bed", "zero.bed", 0.0),
])
def test_metric_empty_file(test_data, capfd, jaccard, a, b, value):
    a_path = test_data("metrics/" + a)
    b_path = test_data("metrics/" + b)

    if jaccard:
        metric = _run_metric_jaccard
    else:
        metric = _run_metric_intersection

    args = [1, 3, "foo"]
    assert (value, *args) == metric(a_path, b_path, *args)

    out, err = capfd.readouterr()
    if not jaccard and a == "zero.bed":
        error = "Warning: Bed file is empty: {}\n".format(
            test_data("metrics/zero.bed")
        )
    else:
        error = ""
    assert error == err


@pytest.mark.parametrize("name,color", [
    ("od", "b"),
    ("od1", "b"),
    ("od_od1", "b"),
    ("ods", "b"),
    ("foo_od2_boo", "b"),
    ("foyd_od", "b"),
    ("OD1", "b"),
    ("fo_OD2", "b"),
    ("fo_OD_YD1", "b"),

    ("yd", "r"),
    ("yd1", "r"),
    ("yd_yd1", "r"),
    ("yds", "r"),
    ("foo_yd2_boo", "r"),
    ("food_yd", "r"),
    ("YD1", "r"),
    ("fo_YD2", "r"),
    ("fo_YD_OD1", "r"),

    ("fo.OD2", "gray"),
    ("fo.YD2", "gray"),
])
def test_color_annotator_age(name, color):
    res = color_annotator_age(name)
    assert 1 == len(res)
    assert "age" == res[0][0]
    assert color == res[0][1]


# To replace images use:
# cd ~/work/washu/test/testdata/metrics
# find . -name "test_plot*.png" | xargs -I fname bash -c "echo fname |
#    sed s/test_plot.*#//g | sed s/.$PLATFORM_SRC$.png//g |
#    xargs -I nname cp fname nname.$PLATFORM_TARGET$.png"
def assert_image(expected_path, actual_path):
    os_platform = platform.system().lower()
    assert os_platform in ["linux", "darwin"], \
        "Unsupported platform: " + os_platform

    if os_platform == "linux":
        os_specific_path = Path(expected_path)
    else:
        os_specific_path \
            = Path(expected_path).with_suffix(".{}.png".format(os_platform))

    try:
        with open(str(os_specific_path), "rb") as ef:
            with open(actual_path, "rb") as af:
                expected = ef.read()
                actual = af.read()
                assert len(expected) == len(actual)
                assert expected == actual

    except (AssertionError, FileNotFoundError):
        prefix = Path(actual_path).parent.name
        exp_name = os_specific_path.name

        if "TC_CHECKOUT_DIR" in os.environ:
            tc_checkout_dir = os.environ['TC_CHECKOUT_DIR']
            export_dir = Path(tc_checkout_dir) / "testsArtifacts"
        else:
            export_dir = os_specific_path.parent

        export_dir.mkdir(exist_ok=True, parents=True)
        export_path = str(export_dir / (prefix + "#" + exp_name))
        shutil.copy(actual_path, export_path)
        print("Mismatched image copied to:", str(export_path),
              file=sys.stderr)
        raise


def test_plot_metric_heatmap(tmp_dir, test_data):
    df = pd.DataFrame.from_csv(test_data("metrics/metric1.csv"))

    result = tmp_dir + "foo.png"
    plot_metric_heatmap("My title", df, save_to=result)

    assert_image(test_data("metrics/img1.png"), result)


@pytest.mark.parametrize("fname,col,row", [
    ("img1.png", False, False),
    ("img1_l-c.png", True, False),
    ("img1_l-r.png", False, True),
    ("img1_l-cr.png", True, True),
])
def test_plot_metric_heatmap_label_fun(tmp_dir, test_data, fname, col, row):
    df = pd.DataFrame.from_csv(test_data("metrics/metric1.csv"))

    expected = test_data("metrics/" + fname)
    result = tmp_dir + "/foo.png"

    col_label_fun = None if not col else lambda x: "[c]" + x
    row_label_fun = None if not row else lambda x: "[r]" + x
    plot_metric_heatmap("My title", df, save_to=result,
                        col_label_converter=col_label_fun,
                        row_label_converter=row_label_fun)

    assert_image(expected, result)


@pytest.mark.parametrize("fname,adjustments", [
    ("img1.png", None),
    ("img1.png", dict()),
    ("img1l.png", dict(left=0.5)),
    ("img1r.png", dict(right=0.5)),
    ("img1t.png", dict(top=0.5)),
    ("img1b.png", dict(bottom=0.5)),
])
def test_plot_metric_heatmap_label_fun(tmp_dir, test_data, fname, adjustments):
    df = pd.DataFrame.from_csv(test_data("metrics/metric1.csv"))

    expected = test_data("metrics/" + fname)
    result = tmp_dir + "/foo.png"

    kw = {}
    if adjustments is not None:
        kw = {**kw, 'adjustments': adjustments}

    plot_metric_heatmap("My title", df, save_to=result, **kw)

    assert_image(expected, result)


@pytest.mark.parametrize("fdf,fname,col,row", [
    ("metric1.csv", "img1_c.png", False, False),
    ("metric1.csv", "img1_c-c.png", True, False),
    ("metric1.csv", "img1_c-r.png", False, True),
    ("metric1.csv", "img1_c-rl.png", True, True),
    ("empty.csv", "empty.png", False, False),
    ("empty.csv", "empty.png", True, False),
    ("empty.csv", "empty.png", False, True),
    ("empty.csv", "empty.png", True, True),
])
def test_plot_metric_heatmap_col_fun(tmp_dir, test_data, fdf, fname, col, row):
    df = pd.DataFrame.from_csv(test_data("metrics/" + fdf))
    if len(df) > 0:
        df.columns = ["ODS", "boo", "YDS", "YD20"]
        df.index = ["OD1", "foo", "boo_YD1"]

    expected = test_data("metrics/" + fname)
    result = "{}/{}".format(tmp_dir, fname)

    col_col_fun = None if not col else color_annotator_age
    row_col_fun = None if not row else color_annotator_age
    plot_metric_heatmap("My title", df, save_to=result,
                        col_color_annotator=col_col_fun,
                        row_color_annotator=row_col_fun)

    assert_image(expected, result)


@pytest.mark.parametrize("fdf,fname,col,row", [
    ("metric1.csv", "img_multi_c.png", False, False),
    ("metric1.csv", "img_multi_c-c.png", True, False),
    ("metric1.csv", "img_multi_c-r.png", False, True),
    ("metric1.csv", "img_multi_c_c-r.png", True, True),
])
def test_plot_multiple_col_fun(tmp_dir, test_data, fdf, fname, col, row):
    df = pd.DataFrame.from_csv(test_data("metrics/" + fdf))

    expected = test_data("metrics/" + fname)
    result = "{}/{}".format(tmp_dir, fname)

    def col_fun(s):
        def inner(label):
            has_a_color = "green" if ("a" in label) else "red"
            has_str = "white" if (s in label) else "gray"
            return (("has_a", has_a_color), ("has_" + s, has_str))
        return inner

    col_col_fun = None if not col else col_fun("1")
    row_col_fun = None if not row else col_fun("2")

    plot_metric_heatmap("My title 2", df, save_to=result,
                        col_color_annotator=col_col_fun,
                        row_color_annotator=row_col_fun)

    assert_image(expected, result)


@pytest.mark.parametrize("label,value", [
    ("OD1", "OD1"),
    ("OD1_foo.broadPeak", "OD1_macs2"),
    ("OD1_foo.narrowPeak", "OD1_macs2"),
    ("OD1_foo-island.bed", "OD1_sicer"),

    ("ODS_foo.broadPeak", "ODS_macs2"),
    ("boo_ODS_foo.broadPeak", "ODS_macs2"),
    ("boo_OD12_foo.broadPeak", "OD12_macs2"),

    ("ODS_foo_peaks.bed", "ODS_zinbra"),
    ("boo_ODS_foo_peaks.bed", "ODS_zinbra"),
    ("boo_OD12_foo_peaks.bed", "OD12_zinbra"),

    ("YDS_boo_consensus_zinbra.bed", "YDS_consensus_zinbra"),
    ("boo_zinbra_YDS_boo_consensus.bed", "YDS_consensus_zinbra"),
    ("boo_zinbra_boo_consensus.bed", "consensus_zinbra"),

    ("YDS_boo_consensus_macs2.bed", "YDS_consensus_macs2"),
    ("boo_macs2_YDS_boo_consensus.bed", "YDS_consensus_macs2"),
    ("boo_macs2_boo_consensus.bed", "consensus_macs2"),
    ("not_tool.bed", "not_tool.bed"),
    ("H3K4me1_consensus.bed", "consensus"),
    ("H3K36me3_consensus.bed", "consensus"),
])
def test_label_converter_donor_and_tool(label, value):
    assert value == label_converter_donor_and_tool(label)


@pytest.mark.parametrize("label,value", [
    ("label1", "(('a', '1'), ('b', '0'), ('c', '0'))"),
    ("label13", "(('a', '1'), ('b', '0'), ('c', '3'))"),
    ("label2", "(('a', '0'), ('b', '2'), ('c', '0'))"),
])
def test_color_annotator_chain(label, value):
    ann1 = lambda l: (("a", "1" if "1" in l else "0"),)
    ann2 = lambda l: (("b", "2" if "2" in l else "0"),)
    ann3 = lambda l: (("c", "3" if "3" in l else "0"),)
    assert value == str(color_annotator_chain(ann1, ann2, ann3)(label))


@pytest.mark.parametrize("label,value", [
    ("label", "(('a', '0'),)"),
    ("label1", "(('a', '1'),)"),
])
def test_color_annotator_chain_single(label, value):
    ann1 = lambda l: (("a", "1" if "1" in l else "0"),)
    assert value == str(color_annotator_chain(ann1)(label))


@pytest.mark.parametrize("label,data_type,value", [
    ("foo_OD1_boo", "H3K4me1", "white"),
    ("foo_OD7_boo", "H3K4me1", "lightgray"),
    ("OD11_boo", "H3K4me1", "g"),
    ("foo_OD11boo", "H3K4me1", "white"),
    ("YD1", "H3K27ac", "lightgray"),
    ("OD20", "methylation", "g"),
])
def test_color_annotator_outlier(test_data, label, data_type, value):
    df = pd.read_csv(test_data("metrics/Y20O20.outliers.csv"),
                     delimiter="\t", skiprows=1, index_col="donor")
    res = color_annotator_outlier(df, data_type)(label)
    assert len(res) == 1
    k, v = res[0]
    assert k == "outlier"
    assert value == v
