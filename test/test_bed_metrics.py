import platform, os
import shutil
from pathlib import Path
import pandas as pd
import pytest
from test.fixtures import test_data, tmp_dir

from reports.bed_metrics import bed_metric_table, heatmap_donor_color_fun, \
    plot_metric_heatmap


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
def test_heatmap_donor_color_fun(name, color):
    res = heatmap_donor_color_fun()(name)
    assert 1 == len(res)
    assert "age" == res[0][0]
    assert color == res[0][1]


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

    except AssertionError:
        prefix = Path(actual_path).parent.name
        exp_name = os_specific_path.name

        if "TC_CHECKOUT_DIR" in os.environ:
            tc_checkout_dir = os_platform.environ['TC_CHECKOUT_DIR']
            export_dir = Path(tc_checkout_dir) / "testsArtifacts"
        else:
            export_dir = os_specific_path.parent

        export_dir.mkdir(exist_ok=True, parents=True)
        shutil.copy(actual_path, str(export_dir / (prefix + "-" + exp_name)))
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
                        col_label_fun=col_label_fun,
                        row_label_fun=row_label_fun)

    assert_image(expected, result)


@pytest.mark.parametrize("fname,col,row", [
    ("img1_c.png", False, False),
    ("img1_c-c.png", True, False),
    ("img1_c-r.png", False, True),
    ("img1_c-rl.png", True, True),
])
def test_plot_metric_heatmap_col_fun(tmp_dir, test_data, fname, col, row):
    df = pd.DataFrame.from_csv(test_data("metrics/metric1.csv"))
    df.columns = ["ODS", "boo", "YDS", "YD20"]
    df.index = ["OD1", "foo", "boo_YD1"]

    expected = test_data("metrics/" + fname)
    result = "{}/{}".format(tmp_dir, fname)

    col_col_fun = None if not col else heatmap_donor_color_fun()
    row_col_fun = None if not row else heatmap_donor_color_fun()
    plot_metric_heatmap("My title", df, save_to=result,
                        col_color_fun=col_col_fun,
                        row_color_fun=row_col_fun)

    assert_image(expected, result)
