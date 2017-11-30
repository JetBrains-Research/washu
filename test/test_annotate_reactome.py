import pytest
import os
import subprocess
from test.fixtures import test_data, tmp_path
from reports.annotate_reactome import extract_reactome_id, fetch_title
from pipeline_utils import run, PROJECT_ROOT_PATH

@pytest.mark.parametrize("label,rid", [
    ("foo", None),
    ("2644605", None),
    ("R-HSA-2644605", "2644605"),
    ("r-hsa-2644605", "2644605"),
    ("foo_R-HSA-2644605_boo", "2644605"),
    ("foo_r-hsa-2644605_boo", "2644605"),
    ("foo.R-HSA-2644605-boo", "2644605"),
    ("foo-R-HSA-2644605.boo", "2644605"),
    ("foo-R-HSA-2644605boo", "2644605"),
    ("http://www.reactome.org/content/detail/R-HSA-2644605", "2644605"),
])
def test_extract_reactome_id(label, rid):
    assert rid == extract_reactome_id(label)


@pytest.mark.parametrize("label,title", [
    ("foo", "N/A"),
    ("foo.R-HSA-2644605-boo", "FBXW7 Mutants and NOTCH1 in Cancer"),
])
def test_fetch_title(label, title):
    assert title == fetch_title(label)


def test_cli_help(capfd):
    run("python", "{}/reports/annotate_reactome.py".format(PROJECT_ROOT_PATH),
        "-h")
    output = capfd.readouterr()
    assert """python {}/reports/annotate_reactome.py -h
usage: annotate_reactome.py [-h] [-F FS] [-o PATH] PATH INT_OR_STR

Annotates table containing http://www.reactome.org pathways ids (R-HSA-nnnnn)
with pathways titles

positional arguments:
  PATH        Table path
  INT_OR_STR  Index (>=1) or name of column containing Reactome pathway id. If
              column is and index column w/o name, use ''

optional arguments:
  -h, --help  show this help message and exit
  -F FS       Field separator, auto-detect if not specified (default: None)
  -o PATH     Output path (default: None)
""".format(PROJECT_ROOT_PATH) == output[0]
    assert "" == output[1]


@pytest.mark.parametrize("file,args,result_stdout,result_stderr,fresult", [
    ("empty.csv", ["1"], "", "File is empty: ", None),
    ("no_header_comma.csv", ["1"], "FBXW7 Mutants and NOTCH1 in Cancer", "", None),
    ("no_header_comma.csv", ["0"], "", "Column index should be >= 1, but was: 0", None),
    ("no_header_comma.csv", ["1"], "", "", "no_header_comma.result.txt"),
    ("header_index_comma.csv", ["1"], "FBXW7 Mutants and NOTCH1 in Cancer", "", None),
    ("header_index_comma.csv", ["1"], "", "", "header_index_comma.result1.txt"),
    ("header_index_comma.csv", ["''"], "", "", "header_index_comma.result2.txt"),
    ("header_colname_comma.csv", ["1"], "FBXW7 Mutants and NOTCH1 in Cancer", "", None),
    ("header_colname_comma.csv", ["1"], "", "", "header_colname_comma.result1.txt"),
    ("header_colname_comma.csv", ["data"], "FBXW7 Mutants and NOTCH1 in Cancer", "", None),
    ("header_colname_comma.csv", ["data"], "", "", "header_colname_comma.result2.txt"),
    ("header_colname_comma.csv", ["data", "-F','"], "FBXW7 Mutants and NOTCH1 in Cancer", "",
     None),
    ("header_colname_tab.tsv", ["data"], "", "", "header_colname_comma.result2.txt"),
    ("header_colname_tab.tsv", ["data", "-F'\t'"], "", "", "header_colname_tab.result.txt"),
    ("header_colname_comma.csv", ["data", "-F' '"], "",
     "KeyError: 'the label [data] is not in the [columns]'", None),
])
def test_foo(test_data, tmp_path, capfd, file, args, result_stdout, result_stderr, fresult):
    input_path = str(test_data("reactome/" + file))

    if fresult:
        args.append("-o " + str(tmp_path / "result.txt"))

    run(
        "python", "{}/reports/annotate_reactome.py".format(PROJECT_ROOT_PATH),
        input_path,
        *args
    )
    output = capfd.readouterr()
    print("==== stdout =====")
    print(output[0])
    print("==== stderr =====")
    print(output[1])
    print("=================")
    assert result_stdout in output[0]
    assert result_stderr in output[1]

    if fresult:
        with open(test_data("reactome/" + fresult), 'r') as ef:
            expected = ef.read()
            with open(str(tmp_path / "result.txt"), 'r') as af:
                actual = af.read()
                assert expected == actual
            # Dfs:
    # By Name, by IDX
    # with / w/o header
    # 0   H3K27ac  http://www.reactome.org/content/detail/R-HSA-2...  transcript
    # 1   H3K27ac  http://www.reactome.org/content/detail/R-HSA-2...  transcript