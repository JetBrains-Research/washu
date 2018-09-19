__author__ = 'Roman Chernyatchik'
__email__ = 'roman.chernyatchik@jetbrains.com'

import pytest

from test.fixtures import test_data, tmp_path
from pipeline_utils import run, PROJECT_ROOT_PATH
from scripts.annotate_reactome import extract_reactome_id, fetch_title


@pytest.mark.parametrize("label,rid", [
    ("foo", None),
    ("2644605", None),
    ("R-HSA-2644605", "R-HSA-2644605"),
    ("r-hsa-2644605", "R-HSA-2644605"),
    ("foo_R-HSA-2644605_boo", "R-HSA-2644605"),
    ("foo_r-hsa-2644605_boo", "R-HSA-2644605"),
    ("foo.R-HSA-2644605-boo", "R-HSA-2644605"),
    ("foo-R-HSA-2644605.boo", "R-HSA-2644605"),
    ("foo-R-HSA-2644605boo", "R-HSA-2644605"),
    ("http://www.reactome.org/content/detail/R-HSA-2644605", "R-HSA-2644605"),
])
def test_extract_reactome_id(label, rid):
    assert rid == extract_reactome_id(label)


@pytest.mark.parametrize("label,title,offline", [
    ("foo", "N/A", False),
    ("foo", "N/A", True),
    ("foo.R-HSA-2644605-boo", "FBXW7 Mutants and NOTCH1 in Cancer", False),
    ("foo.R-HSA-2644605-boo", "FBXW7 Mutants and NOTCH1 in Cancer", True),
])
def test_fetch_title(test_data, label, title, offline):
    kw = {}
    if offline:
        kw['pathways_db_path'] = test_data("reactome/ReactomePathways.txt")

    assert title == fetch_title(label, **kw)[0]


@pytest.mark.parametrize("offline", [True, False])
def test_fetch_titles(test_data, offline):
    kw = {}
    if offline:
        kw['pathways_db_path'] = test_data("reactome/ReactomePathways.txt")

    assert ['FBXW7 Mutants and NOTCH1 in Cancer', '2-LTR circle formation'] == fetch_title(
        "R-HSA-2644605", "R-HSA-164843", **kw
    )


def test_cli_help(capfd):
    run("python", "{}/scripts/annotate_reactome.py".format(PROJECT_ROOT_PATH),
        "-h")
    output = capfd.readouterr()
    assert """python {}/scripts/annotate_reactome.py -h
usage: annotate_reactome.py [-h] [--db PATH] [-F FS] [-o PATH] PATH INT_OR_STR

Annotates table containing http://www.reactome.org pathways ids (R-HSA-nnnnn)
with pathways titles

positional arguments:
  PATH        Table path
  INT_OR_STR  Index (>=1) or name of column containing Reactome pathway id. If
              column is and index column w/o name, use ''

optional arguments:
  -h, --help  show this help message and exit
  --db PATH   ReactomePathways.txt file path (default: None)
  -F FS       Field separator, auto-detect if not specified (default: None)
  -o PATH     Output path (default: None)
""".format(PROJECT_ROOT_PATH) == output[0]
    assert "" == output[1]


@pytest.mark.parametrize("file,args,result_stdout,result_stderr,fresult,offline", [
    # offline
    ("empty.csv", ["1"], "", "File is empty: ", None, True),
    ("no_header_comma.csv", ["1"], "FBXW7 Mutants and NOTCH1 in Cancer", "", None, True),
    ("no_header_comma.csv", ["0"], "", "Column index should be >= 1, but was: 0", None, True),
    ("no_header_comma.csv", ["1"], "", "", "no_header_comma.result.txt", True),
    ("header_index_comma.csv", ["1"], "FBXW7 Mutants and NOTCH1 in Cancer", "", None, True),
    ("header_index_comma.csv", ["1"], "", "", "header_index_comma.result1.txt", True),
    ("header_index_comma.csv", ["''"], "", "", "header_index_comma.result2.txt", True),
    ("header_colname_comma.csv", ["1"], "FBXW7 Mutants and NOTCH1 in Cancer", "", None, True),
    ("header_colname_comma.csv", ["1"], "", "", "header_colname_comma.result1.txt", True),
    ("header_colname_comma.csv", ["data"], "FBXW7 Mutants and NOTCH1 in Cancer", "", None, True),
    ("header_colname_comma.csv", ["data"], "", "", "header_colname_comma.result2.txt", True),
    ("header_colname_comma.csv", ["data", "-F','"], "FBXW7 Mutants and NOTCH1 in Cancer", "",
     None, True),
    ("header_colname_tab.tsv", ["data"], "", "", "header_colname_comma.result2.txt", True),
    ("header_colname_tab.tsv", ["data", "-F'\t'"], "", "", "header_colname_tab.result.txt", True),
    ("header_colname_comma.csv", ["data", "-F' '"], "",
     "KeyError: 'the label [data] is not in the [columns]'", None, True),

    # using web access:
    ("header_colname_comma.csv", ["1"], "FBXW7 Mutants and NOTCH1 in Cancer", "", None, False),
])
def test_cli(test_data, tmp_path, capfd, file, args, result_stdout, result_stderr, fresult,
             offline):

    pathways_db_path = test_data("reactome/ReactomePathways.txt") if offline else None
    input_path = str(test_data("reactome/" + file))

    if fresult:
        args.append("-o " + str(tmp_path / "result.txt"))

    if pathways_db_path:
        args.append("--db " + pathways_db_path)

    run(
        "python", "{}/scripts/annotate_reactome.py".format(PROJECT_ROOT_PATH),
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
